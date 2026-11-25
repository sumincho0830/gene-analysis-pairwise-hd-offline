"""Compare three transcriptomics feature-selection strategies through preprocess.

The input is the 18,045-gene, seed-0 replicate-filtered and biologically
gene-filtered matrix from ``experiments/03_gene_filtered``.  Each network is
analysed independently with one of three strategies:

1. CPM detection support: CPM >= 1 in
   min(available rows, max(4, ceil(10%))) of the network rows.
   For cell-integrated networks the rule is evaluated within each cell type
   and a gene is retained when any cell type passes, preserving cell-specific
   expression.
2. Top 9,000 genes by between-sample variance.  Nine thousand is the
   conservative upper edge of the observed 8,000--9,000 visual elbow.
3. Top 5,000 genes by between-sample variance, matching the approximate
   proteomics feature budget as a sensitivity analysis.

No log transform, scaling, or imputation is performed.  PyWGCNA 2.2.1 uses
its default TPMcutoff and goodSamplesGenes criteria.  A cut=inf diagnostic is
run first.  A finite cut is selected only from linkage gaps that are large in
both absolute and relative terms, retain one unique main cluster, remove at
most five percent of rows, and remove only distance-flagged rows.  The final
selected-cut preprocess is then audited and followed by the same ten QC
outputs used elsewhere in the pipeline.

Large matrices are stored as ``.csv.gz`` and the final reusable PyWGCNA object
as ``.p.gz`` because the working drive has limited free space.  Plotly HTML
uses the CDN bundle rather than embedding a multi-megabyte copy per network.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import gc
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform

from format_corrected_transcriptomics_all_fast import (
    SAMPLE_PATTERN,
    normalized_sample_id,
)
from format_files import TRANSCRIPTOMICS_COLUMN_ALIASES, clean_name
from pipeline_provenance import write_run_manifest
from run_transcriptomics_selected_cuts import pywgcna_retain_largest_cluster
from run_transcriptomics_stage1 import (
    PipelineConfig,
    grouped_distance_robust_z,
    pywgcna_pandas_compatibility,
)
from transcriptomics_qc import (
    PROJECT_RE,
    build_source_cohort_map,
    make_sample_metadata,
    read_expression_matrix,
    write_json,
    write_qc_outputs,
)


STRATEGIES = (
    "01_cpm_detection",
    "02_hvg_top9000",
    "03_hvg_top5000",
)
CPM_THRESHOLD = 1.0
CPM_MIN_FRACTION = 0.10
CPM_MIN_SAMPLES = 4
DIAGNOSTIC_ROBUST_Z = 5.0
CUT_MAX_REMOVED_FRACTION = 0.05
CUT_MIN_GAP_PERCENTILE = 0.75
CUT_MIN_RELATIVE_PERCENTILE = 0.75
CUT_MIN_RELATIVE_GAP = 0.05


@dataclass(frozen=True)
class Network:
    threshold: str
    mode: str
    specimen: str
    cell_type: str
    source_path: Path
    output_root: Path

    @property
    def network_id(self) -> str:
        return "_".join(
            (self.threshold, self.mode, self.specimen, self.cell_type)
        )


@dataclass(frozen=True)
class SelectedCut:
    value: float | None
    reason: str
    removed_ids: tuple[str, ...]
    lower_height: float | None = None
    upper_height: float | None = None
    gap: float | None = None
    gap_percentile: float | None = None
    relative_gap: float | None = None
    relative_gap_percentile: float | None = None
    balanced_score: float | None = None

    @property
    def finite(self) -> bool:
        return bool(self.value is not None and np.isfinite(self.value))

    @property
    def serialized(self) -> float | str:
        return float(self.value) if self.finite else "inf"


def safe_write_json(path: Path, payload: dict[str, Any]) -> None:
    write_json(path, payload)


def write_matrix_gz(matrix: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, index=True, index_label="sample_id", compression="gzip")


def matrices_equal(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    return (
        list(left.index.astype(str)) == list(right.index.astype(str))
        and list(left.columns.astype(str)) == list(right.columns.astype(str))
        and np.array_equal(
            left.to_numpy(dtype=float), right.to_numpy(dtype=float)
        )
    )


def link_or_write_matrix_gz(
    matrix: pd.DataFrame,
    path: Path,
    *,
    link_source: Path | None = None,
    link_source_matrix: pd.DataFrame | None = None,
) -> None:
    """Hard-link an identical checkpoint; otherwise write a gzip CSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    can_link = (
        link_source is not None
        and link_source.exists()
        and link_source_matrix is not None
        and matrices_equal(matrix, link_source_matrix)
    )
    if not can_link:
        write_matrix_gz(matrix, path)
        return
    if path.exists():
        path.unlink()
    try:
        os.link(link_source, path)
    except OSError:
        shutil.copyfile(link_source, path)


def load_networks(config: PipelineConfig, strategies: Iterable[str]) -> list[Network]:
    summary_path = (
        config.experiment_root
        / "03_gene_filtered/transcriptomics/qc_network_summary.csv"
    )
    summary = pd.read_csv(summary_path)
    required = {"threshold", "mode", "specimen", "cell_type", "filtered_matrix"}
    missing = required - set(summary.columns)
    if missing:
        raise ValueError(f"Missing network-summary columns: {sorted(missing)}")
    if len(summary) != 56:
        raise ValueError(f"Expected 56 gene-filtered networks, found {len(summary)}")
    networks: list[Network] = []
    for strategy in strategies:
        for row in summary.itertuples(index=False):
            source_path = Path(str(row.filtered_matrix))
            if not source_path.exists():
                raise FileNotFoundError(source_path)
            root = (
                config.experiment_root
                / "04_variance_filtered/transcriptomics"
                / str(row.threshold)
                / strategy
                / str(row.mode)
                / str(row.specimen)
            )
            if str(row.mode) == "cell_separate":
                root /= str(row.cell_type)
            networks.append(
                Network(
                    threshold=str(row.threshold),
                    mode=str(row.mode),
                    specimen=str(row.specimen),
                    cell_type=str(row.cell_type),
                    source_path=source_path,
                    output_root=root,
                )
            )
    return networks


def integrated_metadata(sample_ids: pd.Index) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    for row_id in sample_ids.astype(str):
        if "__CELL__" not in row_id:
            raise ValueError(f"Integrated row lacks __CELL__ suffix: {row_id}")
        source_sample, safe_cell = row_id.split("__CELL__", 1)
        safe_cell = safe_cell.removeprefix("__")
        records.append(
            {
                "sample_id_specimen_cell": row_id,
                "source_sample_id": source_sample,
                "cell_type": safe_cell,
            }
        )
    return pd.DataFrame(records)


def network_metadata(
    network: Network,
    matrix: pd.DataFrame,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
) -> pd.DataFrame:
    if network.mode == "cell_integrated":
        return make_sample_metadata(
            matrix.index,
            clinical=clinical,
            cohort_map=cohort_map,
            integrated_metadata=integrated_metadata(matrix.index),
        )
    return make_sample_metadata(
        matrix.index,
        clinical=clinical,
        cohort_map=cohort_map,
        fixed_cell_type=network.cell_type,
    )


def threshold_source_name(threshold: str) -> str:
    return {
        "threshold_5_cells": "Threshold_5_Cells",
        "threshold_10_cells": "Threshold_10_Cells",
    }[threshold]


def safe_cell_to_source(config: PipelineConfig) -> dict[str, str]:
    values = pd.read_csv(
        config.run_root / "data/metadata/cell_types.csv", header=None
    ).iloc[:, 0].astype(str)
    mapping: dict[str, str] = {}
    for display in values:
        safe = clean_name(display)
        source = TRANSCRIPTOMICS_COLUMN_ALIASES.get(display, display)
        if safe in mapping and mapping[safe] != source:
            raise ValueError(f"Ambiguous safe cell type: {safe}")
        mapping[safe] = source
    return mapping


def original_sid_lookup(config: PipelineConfig, threshold: str) -> dict[tuple[str, str], str]:
    root = (
        config.run_root
        / "data/Transcriptomics_Corrected"
        / threshold_source_name(threshold)
        / "Batch_Corrected_CSVs"
    )
    lookup: dict[tuple[str, str], str] = {}
    for path in sorted(root.glob("*/*.csv")):
        match = SAMPLE_PATTERN.match(path.stem)
        if match is None:
            raise ValueError(f"Unrecognized corrected filename: {path.name}")
        normalized, _ = normalized_sample_id(path)
        project = PROJECT_RE.match(normalized)
        if project is None:
            raise ValueError(f"Cannot parse project ID: {normalized}")
        project_id = project.group(1)
        subject = match.group("base").removeprefix(project_id + "_")
        original_sid = f"{subject}_{match.group('specimen')}"
        key = (path.parent.name, normalized)
        if key in lookup and lookup[key] != original_sid:
            raise ValueError(f"Ambiguous original sid for {key}")
        lookup[key] = original_sid
    if len(lookup) != 135:
        raise ValueError(
            f"Expected 135 corrected sample files for {threshold}, found {len(lookup)}"
        )
    return lookup


def cpm_targets_for_threshold(
    config: PipelineConfig,
    base_networks: list[Network],
    matrices: dict[str, pd.DataFrame],
    metadata: dict[str, pd.DataFrame],
) -> tuple[
    dict[tuple[str, str, str], list[tuple[str, str]]],
    dict[tuple[str, str], int],
]:
    """Map one CPM column key to network/group targets.

    Target group is ``__all__`` for cell-separated networks and the safe cell
    name for cell-integrated networks.
    """
    if not base_networks:
        return {}, {}
    threshold = base_networks[0].threshold
    sid_lookup = original_sid_lookup(config, threshold)
    cell_map = safe_cell_to_source(config)
    targets: dict[tuple[str, str, str], list[tuple[str, str]]] = {}
    expected: dict[tuple[str, str], int] = {}
    for network in base_networks:
        frame = metadata[network.network_id]
        expected_groups: dict[str, int] = {}
        for row_id, row in frame.iterrows():
            source_sample = str(row["source_sample_id"])
            cohort = str(row["cohort"])
            original_sid = sid_lookup.get((cohort, source_sample))
            if original_sid is None:
                raise ValueError(
                    f"No original CPM sid for {(cohort, source_sample)}"
                )
            safe_cell = (
                str(row["cell_type"])
                if network.mode == "cell_integrated"
                else network.cell_type
            )
            source_cell = cell_map.get(safe_cell)
            if source_cell is None:
                raise ValueError(f"No CPM cell-type alias for {safe_cell}")
            group = safe_cell if network.mode == "cell_integrated" else "__all__"
            key = (cohort, original_sid, source_cell)
            targets.setdefault(key, []).append((network.network_id, group))
            expected_groups[group] = expected_groups.get(group, 0) + 1
        for group, count in expected_groups.items():
            expected[(network.network_id, group)] = count
        if len(frame) != matrices[network.network_id].shape[0]:
            raise ValueError(f"Metadata/matrix row mismatch: {network.network_id}")
    return targets, expected


def scan_cpm_detection(
    config: PipelineConfig,
    base_networks: list[Network],
    matrices: dict[str, pd.DataFrame],
    metadata: dict[str, pd.DataFrame],
) -> tuple[
    dict[tuple[str, str], np.ndarray],
    dict[tuple[str, str], int],
    list[dict[str, Any]],
]:
    """Stream the six large CPM workbooks once for one threshold."""
    if not base_networks:
        return {}, {}, []
    threshold = base_networks[0].threshold
    genes = matrices[base_networks[0].network_id].columns.astype(str)
    for network in base_networks[1:]:
        if list(matrices[network.network_id].columns.astype(str)) != list(genes):
            raise ValueError(f"Gene universe differs: {network.network_id}")
    gene_position = {gene: index for index, gene in enumerate(genes)}
    targets, expected = cpm_targets_for_threshold(
        config, base_networks, matrices, metadata
    )
    counts = {
        target: np.zeros(len(genes), dtype=np.int32) for target in expected
    }
    matched_counts = {target: 0 for target in expected}
    audit_rows: list[dict[str, Any]] = []
    cpm_root = (
        config.run_root
        / "data/Transcriptomics_Corrected"
        / threshold_source_name(threshold)
        / "CPM_Matrices"
    )
    workbooks = sorted(cpm_root.glob("CPM_*.xlsx"))
    if len(workbooks) != 6:
        raise ValueError(f"Expected six CPM workbooks in {cpm_root}, found {len(workbooks)}")

    for workbook_path in workbooks:
        cohort = workbook_path.stem.removeprefix("CPM_")
        print(f"[CPM scan] {threshold}: {workbook_path.name}", flush=True)
        workbook = load_workbook(workbook_path, read_only=True, data_only=True)
        sheet = workbook.active
        rows = sheet.iter_rows(values_only=True)
        header = next(rows)
        if not header or str(header[0]) != "gene_symbol":
            raise ValueError(f"Invalid CPM header: {workbook_path}")
        target_positions: dict[tuple[str, str], list[int]] = {}
        matched_columns = 0
        for position, value in enumerate(header[1:]):
            column = str(value)
            if "___" not in column:
                raise ValueError(f"Invalid CPM column: {column}")
            sid, cell = column.split("___", 1)
            for target in targets.get((cohort, sid, cell), []):
                target_positions.setdefault(target, []).append(position)
                matched_counts[target] += 1
                matched_columns += 1
        genes_seen = 0
        genes_used = 0
        for row in rows:
            genes_seen += 1
            gene = str(row[0])
            gene_index = gene_position.get(gene)
            if gene_index is None:
                continue
            genes_used += 1
            values = row[1:]
            for target, positions in target_positions.items():
                detected = 0
                for position in positions:
                    value = values[position]
                    if value is not None and float(value) >= CPM_THRESHOLD:
                        detected += 1
                counts[target][gene_index] += detected
        workbook.close()
        audit_rows.append(
            {
                "threshold": threshold,
                "workbook": str(workbook_path),
                "cohort": cohort,
                "cpm_columns": len(header) - 1,
                "matched_target_assignments": matched_columns,
                "genes_seen": genes_seen,
                "gene_filtered_genes_matched": genes_used,
            }
        )

    mismatched = [
        (target, expected[target], matched_counts.get(target, 0))
        for target in expected
        if expected[target] != matched_counts.get(target, 0)
    ]
    if mismatched:
        raise ValueError(f"CPM/current-row mapping mismatch: {mismatched[:10]}")
    return counts, expected, audit_rows


def save_cpm_cache(
    cache_root: Path,
    threshold: str,
    counts: dict[tuple[str, str], np.ndarray],
    expected: dict[tuple[str, str], int],
    audit_rows: list[dict[str, Any]],
) -> None:
    cache_root.mkdir(parents=True, exist_ok=True)
    ordered = sorted(counts)
    npz_path = cache_root / f"{threshold}_cpm_detection_counts.npz"
    np.savez_compressed(
        npz_path,
        **{f"array_{index:03d}": counts[key] for index, key in enumerate(ordered)},
    )
    safe_write_json(
        cache_root / f"{threshold}_cpm_detection_counts.json",
        {
            "threshold": threshold,
            "CPM_threshold": CPM_THRESHOLD,
            "targets": [
                {
                    "array": f"array_{index:03d}",
                    "network_id": key[0],
                    "group": key[1],
                    "available_rows": int(expected[key]),
                }
                for index, key in enumerate(ordered)
            ],
            "audit_rows": audit_rows,
        },
    )


def load_cpm_cache(
    cache_root: Path,
    threshold: str,
    expected: dict[tuple[str, str], int],
) -> tuple[
    dict[tuple[str, str], np.ndarray],
    dict[tuple[str, str], int],
    list[dict[str, Any]],
] | None:
    npz_path = cache_root / f"{threshold}_cpm_detection_counts.npz"
    json_path = cache_root / f"{threshold}_cpm_detection_counts.json"
    if not npz_path.exists() or not json_path.exists():
        return None
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    cached_expected = {
        (str(row["network_id"]), str(row["group"])): int(row["available_rows"])
        for row in payload.get("targets", [])
    }
    if cached_expected != expected or float(payload.get("CPM_threshold", -1)) != CPM_THRESHOLD:
        return None
    archive = np.load(npz_path)
    counts = {
        (str(row["network_id"]), str(row["group"])): archive[str(row["array"])].copy()
        for row in payload["targets"]
    }
    archive.close()
    return counts, cached_expected, list(payload.get("audit_rows", []))


def feature_metrics(matrix: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    variance = matrix.var(axis=0, ddof=0).astype(float)
    median = matrix.median(axis=0)
    mad = matrix.sub(median, axis=1).abs().median(axis=0).astype(float)
    return variance, mad


def select_features(
    network: Network,
    strategy: str,
    matrix: pd.DataFrame,
    variance: pd.Series,
    mad: pd.Series,
    cpm_counts: dict[tuple[str, str], np.ndarray],
    cpm_expected: dict[tuple[str, str], int],
) -> tuple[pd.DataFrame, pd.Index, dict[str, Any]]:
    table = pd.DataFrame(
        {
            "gene_symbol": matrix.columns.astype(str),
            "stored_order_1based": np.arange(1, matrix.shape[1] + 1),
            "variance": variance.reindex(matrix.columns).to_numpy(),
            "mad": mad.reindex(matrix.columns).to_numpy(),
        }
    )
    if strategy == "01_cpm_detection":
        groups = sorted(
            group
            for network_id, group in cpm_expected
            if network_id == network.network_id
        )
        if not groups:
            raise ValueError(f"No CPM groups for {network.network_id}")
        passes: list[np.ndarray] = []
        fractions: list[np.ndarray] = []
        overall_count = np.zeros(matrix.shape[1], dtype=np.int32)
        overall_n = 0
        for group in groups:
            count = cpm_counts[(network.network_id, group)]
            n_rows = cpm_expected[(network.network_id, group)]
            required = min(
                n_rows,
                max(CPM_MIN_SAMPLES, math.ceil(CPM_MIN_FRACTION * n_rows)),
            )
            safe_group = re.sub(r"[^A-Za-z0-9]+", "_", group).strip("_") or "all"
            table[f"cpm_detected_{safe_group}"] = count
            table[f"cpm_rows_{safe_group}"] = n_rows
            table[f"cpm_required_{safe_group}"] = required
            fraction = count / n_rows
            table[f"cpm_fraction_{safe_group}"] = fraction
            passes.append(count >= required)
            fractions.append(fraction)
            overall_count += count
            overall_n += n_rows
        selected_mask = np.logical_or.reduce(passes)
        table["cpm_detected_rows_total"] = overall_count
        table["cpm_rows_total"] = overall_n
        table["cpm_detection_fraction_total"] = overall_count / overall_n
        table["cpm_detection_fraction_best_group"] = np.maximum.reduce(fractions)
        table["selected"] = selected_mask
        table["selection_rank"] = (
            table["cpm_detection_fraction_best_group"]
            .rank(method="first", ascending=False)
            .astype(int)
        )
        details = {
            "rule": "CPM >= 1 in min(available rows, max(4, ceil(10%))) rows",
            "integrated_rule": "pass in any cell type" if network.mode == "cell_integrated" else "not_applicable",
            "groups": groups,
        }
    else:
        top_n = 9000 if strategy == "02_hvg_top9000" else 5000
        order = np.argsort(-table["variance"].to_numpy(), kind="mergesort")
        selected_mask = np.zeros(len(table), dtype=bool)
        selected_mask[order[: min(top_n, len(order))]] = True
        rank = np.empty(len(table), dtype=int)
        rank[order] = np.arange(1, len(table) + 1)
        table["variance_rank"] = rank
        table["selection_rank"] = rank
        table["selected"] = selected_mask
        details = {
            "rule": f"top {top_n:,} genes by between-sample population variance",
            "top_n": top_n,
            "tie_break": "stored gene order through stable mergesort",
        }
    selected = pd.Index(table.loc[table["selected"], "gene_symbol"].astype(str))
    if len(selected) < 4:
        raise ValueError(f"Too few selected genes for {network.network_id}/{strategy}")
    details.update(
        input_genes=int(matrix.shape[1]),
        selected_genes=int(len(selected)),
        excluded_genes=int(matrix.shape[1] - len(selected)),
    )
    return table, selected, details


def plot_feature_selection(
    table: pd.DataFrame,
    strategy: str,
    output_path: Path,
    title: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor="white")
    if strategy == "01_cpm_detection":
        ranked = table.sort_values(
            ["cpm_detection_fraction_best_group", "variance"],
            ascending=[False, False],
            kind="mergesort",
        ).reset_index(drop=True)
        axes[0].plot(
            np.arange(1, len(ranked) + 1),
            ranked["cpm_detection_fraction_best_group"],
            color="#4C78A8",
            linewidth=1.2,
        )
        axes[0].set_xlabel("Gene rank by best-group CPM detection")
        axes[0].set_ylabel("Best-group fraction with CPM >= 1")
        axes[0].set_title("CPM detection support")
        axes[1].scatter(
            table["cpm_detection_fraction_best_group"],
            np.log10(table["variance"].clip(lower=0) + 1e-12),
            c=np.where(table["selected"], "#2A9D8F", "#B8B8B8"),
            s=7,
            alpha=0.55,
            linewidths=0,
        )
        axes[1].set_xlabel("Best-group fraction with CPM >= 1")
        axes[1].set_ylabel("log10(between-sample variance)")
        axes[1].set_title("Detection support and corrected-value variance")
    else:
        ranked = table.sort_values("variance_rank", kind="mergesort")
        cutoff = int(table["selected"].sum())
        axes[0].plot(
            ranked["variance_rank"],
            ranked["variance"],
            color="#4C78A8",
            linewidth=1.2,
        )
        axes[0].axvline(cutoff, color="#D62728", linestyle="--", linewidth=1.2)
        axes[0].set_xlabel("Variance rank")
        axes[0].set_ylabel("Between-sample variance")
        axes[0].set_title(f"Variance-ranked genes; cutoff={cutoff:,}")
        axes[1].plot(
            ranked["variance_rank"],
            np.log10(ranked["variance"].clip(lower=0) + 1e-12),
            color="#7A5195",
            linewidth=1.2,
        )
        axes[1].axvline(cutoff, color="#D62728", linestyle="--", linewidth=1.2)
        axes[1].set_xlabel("Variance rank")
        axes[1].set_ylabel("log10(between-sample variance)")
        axes[1].set_title("Variance rank on log scale")
    figure.suptitle(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def percentile_rank(values: np.ndarray, value: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(np.mean(finite <= value))


def metadata_guard(metrics: pd.DataFrame, removed_ids: set[str]) -> str:
    warnings: list[str] = []
    for column in ("cohort", "diagnosis", "sex"):
        if column not in metrics:
            continue
        groups = metrics.groupby(column, dropna=False).groups
        exhausted = [
            str(value)
            for value, indices in groups.items()
            if set(map(str, indices)) and set(map(str, indices)) <= removed_ids
        ]
        if exhausted:
            warnings.append(f"removes_all_{column}=" + "|".join(sorted(exhausted)))
    return ";".join(warnings)


def choose_cut(
    linkage: np.ndarray,
    sample_ids: list[str],
    metrics: pd.DataFrame,
) -> tuple[SelectedCut, pd.DataFrame]:
    heights = linkage[:, 2]
    gaps = np.diff(heights)
    upper = heights[1:]
    relative = np.divide(
        gaps,
        upper,
        out=np.zeros_like(gaps, dtype=float),
        where=upper > 0,
    )
    positive_gaps = gaps[gaps > 1e-12]
    positive_relative = relative[gaps > 1e-12]
    flagged = set(
        metrics.index[
            metrics["diagnostic_outlier_candidate"].fillna(False)
        ].astype(str)
    )
    max_removed = max(1, math.floor(len(sample_ids) * CUT_MAX_REMOVED_FRACTION))
    rows: list[dict[str, Any]] = []
    for index, gap in enumerate(gaps):
        if gap <= 1e-12:
            continue
        lower_height = float(heights[index])
        upper_height = float(heights[index + 1])
        cut = (lower_height + upper_height) / 2.0
        labels = fcluster(linkage, t=cut, criterion="distance")
        unique, counts = np.unique(labels, return_counts=True)
        largest = int(counts.max())
        unique_largest = int(np.sum(counts == largest)) == 1
        main_label = int(unique[int(np.argmax(counts))]) if unique_largest else -1
        removed = tuple(
            sorted(
                sample_ids[position]
                for position, label in enumerate(labels)
                if int(label) != main_label
            )
        ) if unique_largest else tuple(sample_ids)
        gap_percentile = percentile_rank(positive_gaps, float(gap))
        relative_gap = float(relative[index])
        relative_percentile = percentile_rank(positive_relative, relative_gap)
        guard = metadata_guard(metrics, set(removed)) if removed else ""
        all_flagged = bool(removed) and set(removed) <= flagged
        eligible = (
            unique_largest
            and 0 < len(removed) <= max_removed
            and all_flagged
            and not guard
            and gap_percentile >= CUT_MIN_GAP_PERCENTILE
            and relative_percentile >= CUT_MIN_RELATIVE_PERCENTILE
            and relative_gap >= CUT_MIN_RELATIVE_GAP
        )
        rows.append(
            {
                "merge_index_lower": index,
                "lower_height": lower_height,
                "upper_height": upper_height,
                "cut_midpoint": cut,
                "absolute_gap": float(gap),
                "absolute_gap_percentile": gap_percentile,
                "relative_gap": relative_gap,
                "relative_gap_percentile": relative_percentile,
                "balanced_score": min(gap_percentile, relative_percentile),
                "unique_largest_cluster": unique_largest,
                "removed_count": len(removed),
                "max_removed_allowed": max_removed,
                "removed_sample_ids": "; ".join(removed),
                "all_removed_are_distance_flagged": all_flagged,
                "metadata_guard": guard,
                "eligible": eligible,
            }
        )
    table = pd.DataFrame(rows)
    eligible = table.loc[table["eligible"]].copy() if len(table) else table
    if eligible.empty:
        return (
            SelectedCut(
                value=None,
                reason=(
                    "no eligible branch gap met absolute, relative, flagged-row, "
                    "maximum-removal, and metadata guards"
                ),
                removed_ids=(),
            ),
            table,
        )
    best = eligible.sort_values(
        [
            "balanced_score",
            "absolute_gap_percentile",
            "relative_gap_percentile",
            "absolute_gap",
            "relative_gap",
            "removed_count",
        ],
        ascending=[False, False, False, False, False, True],
        kind="mergesort",
    ).iloc[0]
    removed_ids = tuple(
        value.strip()
        for value in str(best["removed_sample_ids"]).split(";")
        if value.strip()
    )
    table["selected"] = False
    table.loc[best.name, "selected"] = True
    return (
        SelectedCut(
            value=float(best["cut_midpoint"]),
            reason="highest balanced absolute-gap and relative-gap support among eligible cuts",
            removed_ids=removed_ids,
            lower_height=float(best["lower_height"]),
            upper_height=float(best["upper_height"]),
            gap=float(best["absolute_gap"]),
            gap_percentile=float(best["absolute_gap_percentile"]),
            relative_gap=float(best["relative_gap"]),
            relative_gap_percentile=float(best["relative_gap_percentile"]),
            balanced_score=float(best["balanced_score"]),
        ),
        table,
    )


def plot_dendrogram_with_cut(
    linkage: np.ndarray,
    labels: list[str],
    selected: SelectedCut,
    output_path: Path,
    title: str,
) -> list[int]:
    width = max(18.0, min(50.0, len(labels) / 8.0))
    figure, axis = plt.subplots(figsize=(width, 9), facecolor="white")
    result = dendrogram(
        linkage,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=5,
        color_threshold=float(selected.value) if selected.finite else 0,
        above_threshold_color="#4C78A8",
        ax=axis,
    )
    if selected.finite:
        axis.axhline(
            float(selected.value), color="#D62728", linewidth=1.3, linestyle="--"
        )
        suffix = f"selected cut={float(selected.value):.6g}"
    else:
        suffix = "selected cut=inf"
    axis.set_title(f"{title} ({suffix})")
    axis.set_xlabel("Samples")
    axis.set_ylabel("Euclidean distance; average linkage")
    figure.tight_layout()
    figure.savefig(output_path, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [int(value) for value in result["leaves"]]


def audit_pre_cut(matrix: pd.DataFrame, model: Any) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    from PyWGCNA import WGCNA

    tpm_keep = (matrix > model.TPMcutoff).any(axis=0).to_numpy()
    after_tpm = matrix.loc[:, tpm_keep]
    with pywgcna_pandas_compatibility():
        good_genes, good_samples, _ = WGCNA.goodSamplesGenes(after_tpm.T)
    good_genes = np.asarray(good_genes, dtype=bool)
    good_samples = np.asarray(good_samples, dtype=bool)
    expected = after_tpm.loc[good_samples, good_genes].copy()
    return expected, tpm_keep, good_genes, good_samples


def assert_matrix_equal(actual: pd.DataFrame, expected: pd.DataFrame, label: str) -> None:
    actual.index = actual.index.astype(str)
    actual.columns = actual.columns.astype(str)
    if list(actual.index) != list(expected.index) or list(actual.columns) != list(expected.columns):
        raise ValueError(f"{label}: IDs differ from audited matrix")
    if not np.array_equal(actual.to_numpy(dtype=float), expected.to_numpy(dtype=float)):
        raise ValueError(f"{label}: values differ from audited matrix")


def save_model_compressed(model: Any, output_root: Path) -> Path:
    output_root.mkdir(parents=True, exist_ok=True)
    model.outputPath = output_root.resolve().as_posix() + "/"
    model.save = True
    model.saveWGCNA()
    raw_path = output_root / f"{model.name}.p"
    # PyWGCNA can return immediately before Windows exposes the newly written
    # pickle to a second process.  Wait briefly before handing it to zstd.
    for _ in range(300):
        if raw_path.exists() and raw_path.stat().st_size > 0:
            break
        time.sleep(0.1)
    if not raw_path.exists():
        raise FileNotFoundError(raw_path)
    zstd = shutil.which("zstd")
    if zstd:
        compressed = raw_path.with_suffix(".p.zst")
        # Feed bytes through standard streams so the MSYS zstd executable
        # never has to resolve this very long Windows path itself.
        with raw_path.open("rb") as source, compressed.open("wb") as target:
            subprocess.run(
                [zstd, "-T0", "-10", "-q", "-c"],
                stdin=source,
                stdout=target,
                check=True,
            )
    else:
        compressed = raw_path.with_suffix(".p.gz")
        with raw_path.open("rb") as source, gzip.open(
            compressed, "wb", compresslevel=6
        ) as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    raw_path.unlink()
    return compressed


def preprocess_network(
    config: PipelineConfig,
    network: Network,
    strategy: str,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    selected_matrix_path: Path,
) -> tuple[pd.DataFrame, dict[str, Any], list[Path]]:
    from PyWGCNA import WGCNA

    preprocess_root = network.output_root / "03_preprocess"
    diagnostic_root = preprocess_root / "01_cut_diagnostic"
    selected_root = preprocess_root / "02_selected_cut"
    diagnostic_root.mkdir(parents=True, exist_ok=True)
    selected_root.mkdir(parents=True, exist_ok=True)
    if matrix.shape[0] < 4:
        summary = {
            "status": "not_run_insufficient_samples",
            "reason": "PyWGCNA default goodSamplesGenes requires at least four samples",
            "input_samples": int(matrix.shape[0]),
            "input_genes": int(matrix.shape[1]),
        }
        safe_write_json(diagnostic_root / "preprocess_summary.json", summary)
        safe_write_json(selected_root / "preprocess_summary.json", summary)
        return matrix.copy(), summary, [
            diagnostic_root / "preprocess_summary.json",
            selected_root / "preprocess_summary.json",
        ]

    sample_info = metadata.loc[matrix.index, [
        "project_id", "patient_id", "cohort", "specimen", "cell_type",
        "diagnosis", "sex", "age",
    ]].copy()
    name = f"TRANSCRIPTOMICS_{network.network_id}_{strategy}_cut_diagnostic"
    diagnostic_model = WGCNA(
        name=name,
        geneExp=matrix,
        sampleInfo=sample_info,
        species="Human",
        save=False,
        outputPath=diagnostic_root.resolve().as_posix() + "/",
        figureType="png",
        cut=float("inf"),
    )
    expected, tpm_keep, good_genes, good_samples = audit_pre_cut(
        matrix, diagnostic_model
    )
    with pywgcna_pandas_compatibility():
        diagnostic_model.preprocess(show=False)
    before_cut = diagnostic_model.datExpr.to_df().copy()
    assert_matrix_equal(before_cut, expected, f"{network.network_id}/{strategy}/cut_inf")
    before_path = diagnostic_root / "matrix_before_cut.csv.gz"
    link_or_write_matrix_gz(
        before_cut,
        before_path,
        link_source=selected_matrix_path,
        link_source_matrix=matrix,
    )

    removed_feature_rows = [
        {
            "gene_symbol": gene,
            "removal_stage": "TPMcutoff",
            "reason": f"no sample value greater than TPMcutoff={diagnostic_model.TPMcutoff}",
        }
        for gene in matrix.columns[~tpm_keep]
    ]
    after_tpm = matrix.loc[:, tpm_keep]
    removed_feature_rows.extend(
        {
            "gene_symbol": gene,
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        }
        for gene in after_tpm.columns[~good_genes]
    )
    removed_features = pd.DataFrame(
        removed_feature_rows,
        columns=["gene_symbol", "removal_stage", "reason"],
    )
    removed_features_path = diagnostic_root / "removed_features.csv"
    removed_features.to_csv(removed_features_path, index=False)
    removed_samples_good = pd.DataFrame(
        {
            "sample_id": after_tpm.index[~good_samples].astype(str),
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        }
    )
    removed_samples_good_path = diagnostic_root / "removed_samples.csv"
    removed_samples_good.to_csv(removed_samples_good_path, index=False)

    condensed = pdist(before_cut.to_numpy(dtype=float), metric="euclidean")
    linkage = WGCNA.hclust(condensed, method="average")
    linkage_path = diagnostic_root / "sample_linkage.csv"
    pd.DataFrame(
        linkage,
        columns=["left_node", "right_node", "distance", "cluster_size"],
    ).rename_axis("merge_index").reset_index().to_csv(linkage_path, index=False)
    distances = squareform(condensed)
    np.fill_diagonal(distances, np.nan)
    metrics = metadata.loc[before_cut.index].copy()
    metrics["mean_pairwise_euclidean_distance"] = np.nanmean(distances, axis=1)
    metrics["median_pairwise_euclidean_distance"] = np.nanmedian(distances, axis=1)
    metrics["max_pairwise_euclidean_distance"] = np.nanmax(distances, axis=1)
    metrics["nearest_neighbor_euclidean_distance"] = np.nanmin(distances, axis=1)
    metrics["median_distance_group_robust_z"] = grouped_distance_robust_z(
        metrics, network.mode == "cell_integrated"
    )
    metrics["diagnostic_distance_threshold"] = DIAGNOSTIC_ROBUST_Z
    metrics["diagnostic_outlier_candidate"] = metrics[
        "median_distance_group_robust_z"
    ].ge(DIAGNOSTIC_ROBUST_Z)
    metrics["candidate_action"] = np.where(
        metrics["diagnostic_outlier_candidate"], "eligible_for_cut_review", "none"
    )
    metrics_path = diagnostic_root / "sample_distance_metrics.csv"
    metrics.reset_index(drop=True).to_csv(metrics_path, index=False)
    flagged_path = diagnostic_root / "outlier_candidates_before_cut.csv"
    metrics.loc[metrics["diagnostic_outlier_candidate"]].reset_index(drop=True).to_csv(
        flagged_path, index=False
    )

    selected_cut, candidates = choose_cut(
        linkage, before_cut.index.astype(str).tolist(), metrics
    )
    candidates_path = diagnostic_root / "cut_candidates.csv.gz"
    candidates.to_csv(candidates_path, index=False, compression="gzip")
    selected_cut_path = diagnostic_root / "selected_cut.json"
    safe_write_json(
        selected_cut_path,
        {
            "selected_cut": selected_cut.serialized,
            "finite_cut_applied": selected_cut.finite,
            "reason": selected_cut.reason,
            "removed_sample_ids": list(selected_cut.removed_ids),
            "lower_height": selected_cut.lower_height,
            "upper_height": selected_cut.upper_height,
            "absolute_gap": selected_cut.gap,
            "absolute_gap_percentile": selected_cut.gap_percentile,
            "relative_gap": selected_cut.relative_gap,
            "relative_gap_percentile": selected_cut.relative_gap_percentile,
            "balanced_score": selected_cut.balanced_score,
            "rules": {
                "max_removed_fraction": CUT_MAX_REMOVED_FRACTION,
                "minimum_absolute_gap_percentile": CUT_MIN_GAP_PERCENTILE,
                "minimum_relative_gap_percentile": CUT_MIN_RELATIVE_PERCENTILE,
                "minimum_relative_gap": CUT_MIN_RELATIVE_GAP,
                "all_removed_must_be_distance_flagged": True,
                "retain_unique_largest_cluster": True,
            },
        },
    )
    dendrogram_path = diagnostic_root / "dendrogram_cut_diagnostic.png"
    leaf_order = plot_dendrogram_with_cut(
        linkage,
        before_cut.index.astype(str).tolist(),
        selected_cut,
        dendrogram_path,
        f"Transcriptomics sample clustering: {network.network_id} / {strategy}",
    )
    leaf_path = diagnostic_root / "dendrogram_leaf_order.csv"
    pd.DataFrame(
        {
            "leaf_position_1based": np.arange(1, len(leaf_order) + 1),
            "sample_id": before_cut.index[leaf_order].astype(str),
        }
    ).to_csv(leaf_path, index=False)
    diagnostic_summary_path = diagnostic_root / "preprocess_summary.json"
    diagnostic_summary = {
        "status": "complete_cut_diagnostic",
        "network_id": network.network_id,
        "strategy": strategy,
        "input_samples": int(matrix.shape[0]),
        "input_genes": int(matrix.shape[1]),
        "TPMcutoff": float(diagnostic_model.TPMcutoff),
        "genes_removed_by_TPMcutoff": int((~tpm_keep).sum()),
        "genes_removed_by_goodSamplesGenes": int((~good_genes).sum()),
        "samples_removed_by_goodSamplesGenes": int((~good_samples).sum()),
        "matrix_before_cut_samples": int(before_cut.shape[0]),
        "matrix_before_cut_genes": int(before_cut.shape[1]),
        "diagnostic_outlier_candidates": int(
            metrics["diagnostic_outlier_candidate"].sum()
        ),
        "automatic_selected_cut": selected_cut.serialized,
        "diagnostic_model_object_saved": False,
    }
    safe_write_json(diagnostic_summary_path, diagnostic_summary)

    if selected_cut.finite:
        final_name = f"TRANSCRIPTOMICS_{network.network_id}_{strategy}_selected_cut"
        final_model = WGCNA(
            name=final_name,
            geneExp=matrix,
            sampleInfo=sample_info,
            species="Human",
            save=True,
            outputPath=selected_root.resolve().as_posix() + "/",
            figureType="png",
            cut=float(selected_cut.value),
        )
        selected_expected_pre, _, _, _ = audit_pre_cut(matrix, final_model)
        assert_matrix_equal(
            selected_expected_pre,
            before_cut.copy(),
            f"{network.network_id}/{strategy}/selected_pre_cut",
        )
        labels = fcluster(linkage, t=float(selected_cut.value), criterion="distance")
        unique, cluster_counts = np.unique(labels, return_counts=True)
        largest = int(cluster_counts.max())
        if int(np.sum(cluster_counts == largest)) != 1:
            raise ValueError("Selected cut does not have a unique largest cluster")
        main_label = int(unique[int(np.argmax(cluster_counts))])
        keep_mask = labels == main_label
        expected_after = before_cut.iloc[keep_mask].copy()
        with pywgcna_pandas_compatibility(), pywgcna_retain_largest_cluster():
            final_model.preprocess(show=False)
        final_matrix = final_model.datExpr.to_df().copy()
        assert_matrix_equal(
            final_matrix,
            expected_after,
            f"{network.network_id}/{strategy}/selected_cut",
        )
    else:
        final_model = diagnostic_model
        final_model.name = f"TRANSCRIPTOMICS_{network.network_id}_{strategy}_selected_cut_inf"
        final_matrix = before_cut.copy()
        keep_mask = np.ones(len(before_cut), dtype=bool)
        labels = np.zeros(len(before_cut), dtype=int)
        cluster_counts = np.asarray([len(before_cut)])

    final_model.analysis_provenance = {
        "stage": "04_variance_filtered_selected_cut",
        "network_id": network.network_id,
        "feature_strategy": strategy,
        "transform": "none",
        "imputation": "none",
        "TPMcutoff": final_model.TPMcutoff,
        "goodSamplesGenes": "PyWGCNA 2.2.1 defaults with pandas compatibility fix",
        "cut": selected_cut.serialized,
        "cut_selection": selected_cut.reason,
        "findModules_run": False,
    }
    object_path = save_model_compressed(final_model, selected_root)
    final_matrix_path = selected_root / "matrix_after_preprocess.csv.gz"
    link_or_write_matrix_gz(
        final_matrix,
        final_matrix_path,
        link_source=before_path,
        link_source_matrix=before_cut,
    )
    selected_dendrogram_path = selected_root / "dendrogram_selected_cut.png"
    plot_dendrogram_with_cut(
        linkage,
        before_cut.index.astype(str).tolist(),
        selected_cut,
        selected_dendrogram_path,
        f"Selected-cut transcriptomics clustering: {network.network_id} / {strategy}",
    )
    assignments = metadata.loc[before_cut.index].copy()
    assignments["cluster_label"] = labels
    assignments["selected_cut"] = (
        float(selected_cut.value) if selected_cut.finite else np.nan
    )
    assignments["kept_by_selected_cut"] = keep_mask
    assignments["cut_action"] = np.where(
        keep_mask, "retained", "removed_by_selected_cut"
    )
    assignments_path = selected_root / "sample_cluster_assignments.csv"
    assignments.reset_index(drop=True).to_csv(assignments_path, index=False)
    removed_by_cut = assignments.loc[~assignments["kept_by_selected_cut"]].copy()
    removed_by_cut_path = selected_root / "removed_samples.csv"
    removed_by_cut.reset_index(drop=True).to_csv(removed_by_cut_path, index=False)
    removed_features_selected_path = selected_root / "removed_features.csv"
    removed_features.to_csv(removed_features_selected_path, index=False)
    selected_summary = {
        "status": "complete_selected_cut",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "network_id": network.network_id,
        "strategy": strategy,
        "input_samples": int(matrix.shape[0]),
        "input_genes": int(matrix.shape[1]),
        "TPMcutoff": float(final_model.TPMcutoff),
        "selected_cut": selected_cut.serialized,
        "cut_selection_reason": selected_cut.reason,
        "finite_cut_applied": selected_cut.finite,
        "distance_flagged_before_cut": int(
            metrics["diagnostic_outlier_candidate"].sum()
        ),
        "samples_removed_by_goodSamplesGenes": int((~good_samples).sum()),
        "samples_removed_by_dendrogram_cut": int((~keep_mask).sum()),
        "samples_after_preprocess": int(final_matrix.shape[0]),
        "genes_after_preprocess": int(final_matrix.shape[1]),
        "model_object": str(object_path),
        "model_object_compression": (
            "zstd; decompress to .p before PyWGCNA.loadWGCNA"
            if object_path.suffix == ".zst"
            else "gzip; decompress to .p before PyWGCNA.loadWGCNA"
        ),
        "findModules_run": False,
    }
    selected_summary_path = selected_root / "preprocess_summary.json"
    safe_write_json(selected_summary_path, selected_summary)
    outputs = [
        before_path,
        removed_features_path,
        removed_samples_good_path,
        linkage_path,
        metrics_path,
        flagged_path,
        candidates_path,
        selected_cut_path,
        dendrogram_path,
        leaf_path,
        diagnostic_summary_path,
        object_path,
        final_matrix_path,
        selected_dendrogram_path,
        assignments_path,
        removed_by_cut_path,
        removed_features_selected_path,
        selected_summary_path,
    ]
    del diagnostic_model, final_model
    gc.collect()
    return final_matrix, selected_summary, outputs


def write_feature_outputs(
    network: Network,
    strategy: str,
    matrix: pd.DataFrame,
    table: pd.DataFrame,
    selected: pd.Index,
    details: dict[str, Any],
) -> tuple[pd.DataFrame, list[Path]]:
    root = network.output_root / "01_feature_selection"
    root.mkdir(parents=True, exist_ok=True)
    table_path = root / "feature_selection.csv.gz"
    table.to_csv(table_path, index=False, compression="gzip")
    selected_path = root / "selected_features.csv"
    excluded_path = root / "excluded_features.csv"
    table.loc[table["selected"]].to_csv(selected_path, index=False)
    table.loc[~table["selected"]].to_csv(excluded_path, index=False)
    figure_path = root / "feature_selection.png"
    plot_feature_selection(
        table,
        strategy,
        figure_path,
        f"Feature selection: {network.network_id} / {strategy}",
    )
    filtered = matrix.loc[:, selected].copy()
    matrix_path = root / "selected_expression_matrix.csv.gz"
    write_matrix_gz(filtered, matrix_path)
    summary_path = root / "feature_selection_summary.json"
    safe_write_json(
        summary_path,
        {
            "network_id": network.network_id,
            "strategy": strategy,
            "source_matrix": str(network.source_path),
            "expression_values": "stored batch-corrected values; no transform, scaling, or imputation",
            **details,
        },
    )
    return filtered, [
        table_path,
        selected_path,
        excluded_path,
        figure_path,
        matrix_path,
        summary_path,
    ]


def plot_global_comparison(summary: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(17, 6), facecolor="white")
    strategies = list(STRATEGIES)
    data = [
        summary.loc[summary["strategy"].eq(strategy), "selected_genes"].dropna()
        for strategy in strategies
    ]
    axes[0].boxplot(data, tick_labels=strategies, showfliers=True)
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_ylabel("Genes selected before PyWGCNA")
    axes[0].set_title("Feature budget")
    data = [
        summary.loc[
            summary["strategy"].eq(strategy), "genes_after_preprocess"
        ].dropna()
        for strategy in strategies
    ]
    axes[1].boxplot(data, tick_labels=strategies, showfliers=True)
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].set_ylabel("Genes after PyWGCNA preprocess")
    axes[1].set_title("Preprocess retention")
    removal = (
        summary.groupby("strategy")["samples_removed_by_dendrogram_cut"]
        .sum()
        .reindex(strategies)
        .fillna(0)
    )
    axes[2].bar(strategies, removal, color=["#4C78A8", "#F28E2B", "#59A14F"])
    axes[2].tick_params(axis="x", rotation=25)
    axes[2].set_ylabel("Total rows removed")
    axes[2].set_title("Automatic dendrogram cuts")
    figure.suptitle("Transcriptomics feature-selection and preprocess comparison")
    figure.tight_layout()
    figure.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def run(
    config: PipelineConfig,
    *,
    strategies: tuple[str, ...],
    thresholds: tuple[str, ...],
    only_network: str | None = None,
    resume: bool = True,
) -> None:
    all_networks = load_networks(config, strategies)
    selected_networks = [
        network
        for network in all_networks
        if network.threshold in thresholds
        and (only_network is None or network.network_id == only_network)
    ]
    if not selected_networks:
        raise ValueError("No networks matched the requested filters")

    clinical = pd.read_csv(config.clinical_path)
    cohort_map = build_source_cohort_map(config.run_root)
    base_by_threshold: dict[str, list[Network]] = {}
    for threshold in thresholds:
        seen: set[str] = set()
        base: list[Network] = []
        for network in selected_networks:
            if network.threshold == threshold and network.network_id not in seen:
                seen.add(network.network_id)
                base.append(network)
        base_by_threshold[threshold] = base

    matrices: dict[str, pd.DataFrame] = {}
    metadata: dict[str, pd.DataFrame] = {}
    for networks in base_by_threshold.values():
        for network in networks:
            matrix = read_expression_matrix(network.source_path)
            matrices[network.network_id] = matrix
            metadata[network.network_id] = network_metadata(
                network, matrix, clinical, cohort_map
            )

    stage_root = config.experiment_root / "04_variance_filtered/transcriptomics"
    stage_root.mkdir(parents=True, exist_ok=True)
    cpm_cache_root = stage_root / "00_shared/cpm_detection_cache"
    cpm_counts: dict[tuple[str, str], np.ndarray] = {}
    cpm_expected: dict[tuple[str, str], int] = {}
    cpm_audit_rows: list[dict[str, Any]] = []
    if "01_cpm_detection" in strategies:
        for threshold, networks in base_by_threshold.items():
            _, expected_contract = cpm_targets_for_threshold(
                config, networks, matrices, metadata
            )
            cached = load_cpm_cache(
                cpm_cache_root, threshold, expected_contract
            )
            if cached is None:
                counts, expected, audit = scan_cpm_detection(
                    config, networks, matrices, metadata
                )
                save_cpm_cache(
                    cpm_cache_root, threshold, counts, expected, audit
                )
            else:
                print(f"[CPM cache] {threshold}", flush=True)
                counts, expected, audit = cached
            cpm_counts.update(counts)
            cpm_expected.update(expected)
            cpm_audit_rows.extend(audit)

    cpm_audit_path = stage_root / "cpm_source_audit.csv"
    pd.DataFrame(cpm_audit_rows).to_csv(cpm_audit_path, index=False)
    rows: list[dict[str, Any]] = []
    for position, network in enumerate(selected_networks, start=1):
        strategy = network.output_root.parts[-4] if network.mode == "cell_separate" else network.output_root.parts[-3]
        # The path-based inference above is intentionally checked against the
        # known strategy names to prevent a silent directory-layout mistake.
        strategy_candidates = [value for value in STRATEGIES if value in network.output_root.parts]
        if len(strategy_candidates) != 1:
            raise ValueError(f"Cannot resolve strategy from {network.output_root}")
        strategy = strategy_candidates[0]
        selected_summary_path = (
            network.output_root / "03_preprocess/02_selected_cut/preprocess_summary.json"
        )
        network_manifest = network.output_root / "run_manifest.json"
        if resume and selected_summary_path.exists() and network_manifest.exists():
            summary = json.loads(selected_summary_path.read_text(encoding="utf-8"))
            feature_summary = json.loads(
                (network.output_root / "01_feature_selection/feature_selection_summary.json").read_text(
                    encoding="utf-8"
                )
            )
            manifest_data = json.loads(network_manifest.read_text(encoding="utf-8"))
            compression = summary.get("model_object_compression")
            if compression:
                storage_note = (
                    "large matrices gzip; final PyWGCNA object "
                    f"{compression.split(';', 1)[0]}; Plotly CDN"
                )
            else:
                storage_note = (
                    "large matrices gzip; final PyWGCNA object not created "
                    "(insufficient samples); Plotly CDN"
                )
            if manifest_data.get("generation_rules", {}).get("storage") != storage_note:
                manifest_data.setdefault("generation_rules", {})["storage"] = storage_note
                safe_write_json(network_manifest, manifest_data)
            rows.append({
                "threshold": network.threshold,
                "strategy": strategy,
                "mode": network.mode,
                "specimen": network.specimen,
                "cell_type": network.cell_type,
                "network_id": network.network_id,
                "status": summary.get("status"),
                "input_samples": summary.get("input_samples"),
                "input_genes": feature_summary.get("input_genes"),
                "selected_genes": feature_summary.get("selected_genes"),
                **{key: summary.get(key) for key in (
                    "genes_after_preprocess", "selected_cut",
                    "distance_flagged_before_cut",
                    "samples_removed_by_dendrogram_cut", "samples_after_preprocess",
                )},
            })
            continue
        print(
            f"[variance-filter {position}/{len(selected_networks)}] "
            f"{strategy} / {network.network_id}",
            flush=True,
        )
        matrix = matrices[network.network_id]
        network_meta = metadata[network.network_id]
        variance, mad = feature_metrics(matrix)
        selection_table, selected, details = select_features(
            network,
            strategy,
            matrix,
            variance,
            mad,
            cpm_counts,
            cpm_expected,
        )
        filtered, feature_paths = write_feature_outputs(
            network,
            strategy,
            matrix,
            selection_table,
            selected,
            details,
        )
        before_root = network.output_root / "02_before_preprocess"
        before_qc = write_qc_outputs(
            filtered,
            network_meta,
            before_root,
            title=f"Variance-filtered transcriptomics: {network.network_id} / {strategy}",
            include_plotlyjs="cdn",
        )
        final_matrix, preprocess_summary, preprocess_paths = preprocess_network(
            config,
            network,
            strategy,
            filtered,
            network_meta,
            network.output_root
            / "01_feature_selection/selected_expression_matrix.csv.gz",
        )
        after_root = network.output_root / "04_after_preprocess"
        after_meta = network_meta.loc[final_matrix.index]
        after_qc = write_qc_outputs(
            final_matrix,
            after_meta,
            after_root,
            title=f"After PyWGCNA preprocess: {network.network_id} / {strategy}",
            include_plotlyjs="cdn",
        )
        write_run_manifest(
            network.output_root / "run_manifest.json",
            pipeline_name="transcriptomics_variance_filter_preprocess_comparison",
            run_root=config.run_root,
            parameters={
                "threshold": network.threshold,
                "strategy": strategy,
                "mode": network.mode,
                "specimen": network.specimen,
                "cell_type": network.cell_type,
                "CPM_detection_threshold": CPM_THRESHOLD if strategy == "01_cpm_detection" else None,
                "CPM_min_fraction": CPM_MIN_FRACTION if strategy == "01_cpm_detection" else None,
                "CPM_min_samples": CPM_MIN_SAMPLES if strategy == "01_cpm_detection" else None,
            },
            generation_rules={
                "input": "03_gene_filtered matrix after CIS exclusion and seed-0 replicate selection",
                "expression_transform": "none",
                "imputation": "none",
                "variance": "population variance ddof=0 within each network",
                "cut_selection": (
                    "midpoint of an eligible linkage gap with both absolute and relative gap percentile >= 0.75, "
                    "relative gap >= 0.05, <=5% removed, all removed distance-flagged, metadata levels retained"
                ),
                "after_preprocess": "ten descriptive/distribution outputs regenerated from final matrix",
                "findModules": "not run",
                "storage": (
                    "large matrices gzip; final PyWGCNA object "
                    f"{preprocess_summary.get('model_object_compression', 'not created (insufficient samples)').split(';', 1)[0]}; "
                    "Plotly CDN"
                ),
            },
            input_paths=[network.source_path, config.clinical_path],
            output_paths=[
                *feature_paths,
                *before_qc.paths,
                *preprocess_paths,
                *after_qc.paths,
            ],
            code_paths=[
                Path(__file__),
                Path(__file__).parent / "transcriptomics_qc.py",
                Path(__file__).parent / "run_transcriptomics_stage1.py",
                Path(__file__).parent / "run_transcriptomics_selected_cuts.py",
            ],
            packages=(
                "numpy", "pandas", "matplotlib", "scipy", "scikit-learn",
                "plotly", "openpyxl", "PyWGCNA", "anndata",
            ),
        )
        rows.append(
            {
                "threshold": network.threshold,
                "strategy": strategy,
                "mode": network.mode,
                "specimen": network.specimen,
                "cell_type": network.cell_type,
                "network_id": network.network_id,
                "status": preprocess_summary.get("status"),
                "input_samples": int(matrix.shape[0]),
                "input_genes": int(matrix.shape[1]),
                "selected_genes": int(filtered.shape[1]),
                "genes_after_preprocess": preprocess_summary.get("genes_after_preprocess"),
                "selected_cut": preprocess_summary.get("selected_cut"),
                "distance_flagged_before_cut": preprocess_summary.get("distance_flagged_before_cut"),
                "samples_removed_by_dendrogram_cut": preprocess_summary.get("samples_removed_by_dendrogram_cut"),
                "samples_after_preprocess": preprocess_summary.get("samples_after_preprocess"),
            }
        )
        gc.collect()

    summary_path = stage_root / "variance_filter_preprocess_summary.csv"
    new_summary = pd.DataFrame(rows)
    if resume and summary_path.exists():
        old = pd.read_csv(summary_path)
        combined = pd.concat([old, new_summary], ignore_index=True)
        combined = combined.drop_duplicates(
            ["threshold", "strategy", "mode", "specimen", "cell_type"],
            keep="last",
        )
    else:
        combined = new_summary
    combined.to_csv(summary_path, index=False)
    comparison_path = stage_root / "variance_filter_preprocess_comparison.png"
    plot_global_comparison(combined, comparison_path)
    write_run_manifest(
        stage_root / "run_manifest.json",
        pipeline_name="all_transcriptomics_variance_filter_preprocess_comparisons",
        run_root=config.run_root,
        parameters={
            "strategies": list(strategies),
            "thresholds": list(thresholds),
            "networks_recorded": int(len(combined)),
            "only_network": only_network,
        },
        generation_rules={
            "strategy_1": "CPM detection support on source CPM, applied to corrected matrix",
            "strategy_2": "top 9,000 corrected-value variance genes",
            "strategy_3": "top 5,000 corrected-value variance genes",
            "preprocess": "cut=inf diagnostic then automatic conservative selected cut",
            "after_preprocess": "descriptive and distribution QC regenerated",
            "stop_gate": "findModules not run",
        },
        input_paths=[
            config.experiment_root
            / "03_gene_filtered/transcriptomics/qc_network_summary.csv",
            config.clinical_path,
        ],
        output_paths=[cpm_audit_path, summary_path, comparison_path],
        code_paths=[Path(__file__)],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "openpyxl"),
    )


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=default_root)
    parser.add_argument("--strategy", action="append", choices=STRATEGIES)
    parser.add_argument(
        "--threshold",
        action="append",
        choices=("threshold_5_cells", "threshold_10_cells"),
    )
    parser.add_argument("--only-network")
    parser.add_argument("--no-resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(
        PipelineConfig(run_root=args.run_root.resolve()),
        strategies=tuple(args.strategy or STRATEGIES),
        thresholds=tuple(
            args.threshold or ("threshold_5_cells", "threshold_10_cells")
        ),
        only_network=args.only_network,
        resume=not args.no_resume,
    )
    print("Completed transcriptomics variance-filter preprocess comparison.")


if __name__ == "__main__":
    main()
