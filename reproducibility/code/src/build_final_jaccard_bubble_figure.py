"""Build final CSF RNA-protein module bubble figures with permutation FDR.

The figure family is intentionally fixed to the analysis requested for the
manuscript candidate:

* specimen: CSF
* transcriptomics mode: cell_separate
* feature strategy: HVG top 9,000
* sample branch: after removal of Ng_MOG1 and Ng_HD3
* cell types: every saved cell type except ``other``
* CPM detection thresholds: 5 cells and 10 cells, reported separately

Dot area is proportional to the saved |kME|-weighted Jaccard similarity.  Dot
colour is the Benjamini-Hochberg FDR q-value from a one-sided permutation test.
For each cell-type network, the RNA gene identities are permuted within that
network's RNA/protein common-gene universe.  This preserves RNA module sizes,
RNA module kME weights, protein module sizes, and protein module kME weights,
while breaking gene-level correspondence between the two omics layers.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

from run_cross_omics_kme_weighted_jaccard import (
    collapse_proteomics_to_gene,
    module_weights,
    normalize_gene,
)


ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = ROOT / "figures" / "final_jaccard_bubble"
STRATEGY = "02_hvg_top9000"
PROTEOMICS_RECORD_ID = "proteomics::CSF::04_ds3_mms30_med0p10"
THRESHOLDS = ("threshold_5_cells", "threshold_10_cells")
CORRECTION_LABELS = {
    "combat_v1_cell_diagnosis": "ComBat v1",
    "combat_v2_cell_diagnosis_tissue": "ComBat v2",
}
DEFAULT_CORRECTION = "combat_v1_cell_diagnosis"
CELL_ORDER = (
    "B_cell",
    "CD4_T_cell",
    "CD8_T_cell",
    "cDCs",
    "Classical_Monocytes",
    "Macrophages",
    "MAIT_gDT_cell",
    "NK_cell",
    "Non_classical_Monocytes",
    "pDCs",
    "Plasmablast",
    "Proliferating_T_cell",
)
CELL_LABELS = {
    "B_cell": "B cell",
    "CD4_T_cell": "CD4 T cell",
    "CD8_T_cell": "CD8 T cell",
    "cDCs": "cDCs",
    "Classical_Monocytes": "Classical monocytes",
    "Macrophages": "Macrophages",
    "MAIT_gDT_cell": "MAIT / γδ T cell",
    "NK_cell": "NK cell",
    "Non_classical_Monocytes": "Non-classical monocytes",
    "pDCs": "pDCs",
    "Plasmablast": "Plasmablast",
    "Proliferating_T_cell": "Proliferating T cell",
}


@dataclass(frozen=True)
class NetworkFiles:
    record_id: str
    threshold: str
    cell_type: str
    rna_kme: Path
    pair_table: Path
    contribution_table: Path
    diagnosis_table: Path | None


def read_dashboard_records(dashboard_root: Path = ROOT) -> list[dict[str, Any]]:
    source = (dashboard_root / "dashboard-data.js").read_text(encoding="utf-8")
    payload = source.split("=", 1)[1].strip().removesuffix(";")
    return json.loads(payload)["records"]


def file_map(record: dict[str, Any], dashboard_root: Path = ROOT) -> dict[str, Path]:
    return {item["id"]: dashboard_root / item["path"] for item in record.get("files", [])}


def select_networks(
    records: list[dict[str, Any]],
    batch_correction: str = DEFAULT_CORRECTION,
    dashboard_root: Path = ROOT,
) -> list[NetworkFiles]:
    selected: list[NetworkFiles] = []
    for record in records:
        dimensions = record.get("dimensions", {})
        cell_type = str(dimensions.get("cell_type", ""))
        if not (
            dimensions.get("threshold") in THRESHOLDS
            and dimensions.get("batch_correction", DEFAULT_CORRECTION) == batch_correction
            and dimensions.get("strategy") == STRATEGY
            and dimensions.get("mode") == "cell_separate"
            and dimensions.get("specimen") == "CSF"
            and cell_type.lower() != "other"
        ):
            continue
        files = file_map(record, dashboard_root)
        required = (
            "after_outlier_official_kme",
            "after_outlier_cross_pair_similarity",
            "after_outlier_cross_all_gene_contributions",
        )
        if not all(key in files for key in required):
            continue
        selected.append(
            NetworkFiles(
                record_id=record["id"],
                threshold=str(dimensions["threshold"]),
                cell_type=cell_type,
                rna_kme=files["after_outlier_official_kme"],
                pair_table=files["after_outlier_cross_pair_similarity"],
                contribution_table=files[
                    "after_outlier_cross_all_gene_contributions"
                ],
                diagnosis_table=files.get(
                    "after_outlier_diagnosis_pairwise_statistics"
                ),
            )
        )
    order = {name: index for index, name in enumerate(CELL_ORDER)}
    return sorted(
        selected,
        key=lambda item: (
            THRESHOLDS.index(item.threshold),
            order.get(item.cell_type, 999),
            item.cell_type,
        ),
    )


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Return BH-FDR adjusted p-values, preserving the original order."""

    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if len(values) == 0:
        return values.copy()
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("p_values must be finite and lie in [0, 1]")
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0, 1)
    return result


def module_arrays(
    genes: list[str],
    modules: dict[str, dict[str, float]],
    module_order: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions = {gene: index for index, gene in enumerate(genes)}
    module_index = np.full(len(genes), -1, dtype=np.int32)
    weight = np.zeros(len(genes), dtype=np.float64)
    sums = np.zeros(len(module_order), dtype=np.float64)
    for index, module in enumerate(module_order):
        for gene, value in modules[module].items():
            position = positions[gene]
            if module_index[position] >= 0:
                raise ValueError(f"Gene {gene} occurs in more than one module")
            module_index[position] = index
            weight[position] = float(value)
            sums[index] += float(value)
    return module_index, weight, sums


def observed_weighted_jaccard(
    rna_module_index: np.ndarray,
    rna_weight: np.ndarray,
    rna_sums: np.ndarray,
    protein_module_index: np.ndarray,
    protein_weight: np.ndarray,
    protein_sums: np.ndarray,
) -> np.ndarray:
    module_count = len(protein_sums)
    mask = (rna_module_index >= 0) & (protein_module_index >= 0)
    flat_index = (
        rna_module_index[mask].astype(np.int64) * module_count
        + protein_module_index[mask].astype(np.int64)
    )
    intersection = np.bincount(
        flat_index,
        weights=np.minimum(rna_weight[mask], protein_weight[mask]),
        minlength=len(rna_sums) * module_count,
    ).reshape(len(rna_sums), module_count)
    union = rna_sums[:, None] + protein_sums[None, :] - intersection
    return np.divide(
        intersection,
        union,
        out=np.zeros_like(intersection),
        where=union > 0,
    )


def permutation_p_values(
    observed: np.ndarray,
    rna_module_index: np.ndarray,
    rna_weight: np.ndarray,
    rna_sums: np.ndarray,
    protein_module_index: np.ndarray,
    protein_weight: np.ndarray,
    protein_sums: np.ndarray,
    *,
    permutations: int,
    seed: int,
) -> np.ndarray:
    """One-sided empirical p-values for unusually large weighted Jaccard."""

    if permutations < 1:
        raise ValueError("permutations must be positive")
    rng = np.random.default_rng(seed)
    exceedances = np.zeros(observed.shape, dtype=np.int64)
    protein_module_count = len(protein_sums)
    for _ in range(permutations):
        permutation = rng.permutation(len(rna_module_index))
        permuted_module = rna_module_index[permutation]
        permuted_weight = rna_weight[permutation]
        mask = (permuted_module >= 0) & (protein_module_index >= 0)
        flat_index = (
            permuted_module[mask].astype(np.int64) * protein_module_count
            + protein_module_index[mask].astype(np.int64)
        )
        intersection = np.bincount(
            flat_index,
            weights=np.minimum(permuted_weight[mask], protein_weight[mask]),
            minlength=len(rna_sums) * protein_module_count,
        ).reshape(observed.shape)
        union = rna_sums[:, None] + protein_sums[None, :] - intersection
        null = np.divide(
            intersection,
            union,
            out=np.zeros_like(intersection),
            where=union > 0,
        )
        exceedances += null >= (observed - 1e-15)
    # The +1 correction prevents zero Monte Carlo p-values.
    return (exceedances + 1.0) / (permutations + 1.0)


def analyse_network(
    network: NetworkFiles,
    protein_collapsed: pd.DataFrame,
    protein_universe: set[str],
    protein_module_order: list[str],
    *,
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    pairs = pd.read_csv(network.pair_table)
    rna = pd.read_csv(network.rna_kme)
    rna["gene_symbol"] = rna["gene_id"].map(normalize_gene)
    common = set(rna["gene_symbol"].dropna().astype(str)) & protein_universe

    protein = protein_collapsed.copy()
    protein["is_raw_unassigned"] = False
    rna_order = pairs["transcript_module"].drop_duplicates().astype(str).tolist()
    pair_protein_order = (
        pairs["proteomics_module"].drop_duplicates().astype(str).tolist()
    )
    if pair_protein_order != protein_module_order:
        raise ValueError(
            f"Protein module order changed for {network.record_id}: "
            f"{pair_protein_order}"
        )
    rna_modules, _ = module_weights(
        rna,
        "gene_symbol",
        "official_assigned_abs_kME",
        common,
        rna_order,
    )
    protein_modules, _ = module_weights(
        protein,
        "gene_symbol",
        "official_assigned_abs_kME",
        common,
        protein_module_order,
    )
    genes = sorted(common)
    rna_index, rna_weight, rna_sums = module_arrays(
        genes, rna_modules, rna_order
    )
    protein_index, protein_weight, protein_sums = module_arrays(
        genes, protein_modules, protein_module_order
    )
    observed = observed_weighted_jaccard(
        rna_index,
        rna_weight,
        rna_sums,
        protein_index,
        protein_weight,
        protein_sums,
    )
    stored = (
        pairs.pivot(
            index="transcript_module",
            columns="proteomics_module",
            values="kme_weighted_jaccard",
        )
        .reindex(index=rna_order, columns=protein_module_order)
        .to_numpy(dtype=float)
    )
    maximum_error = float(np.nanmax(np.abs(observed - stored)))
    if maximum_error > 1e-10:
        raise ValueError(
            f"Stored/recomputed weighted Jaccard mismatch for {network.record_id}: "
            f"max abs error={maximum_error:.3g}"
        )
    p_values = permutation_p_values(
        observed,
        rna_index,
        rna_weight,
        rna_sums,
        protein_index,
        protein_weight,
        protein_sums,
        permutations=permutations,
        seed=seed,
    )
    p_table = pd.DataFrame(
        p_values,
        index=pd.Index(rna_order, name="transcript_module"),
        columns=pd.Index(protein_module_order, name="proteomics_module"),
    ).stack(future_stack=True).rename("permutation_p_value").reset_index()
    result = pairs.merge(
        p_table,
        on=["transcript_module", "proteomics_module"],
        how="left",
        validate="one_to_one",
    )
    result.insert(0, "record_id", network.record_id)
    result.insert(1, "threshold", network.threshold)
    result.insert(2, "cell_type", network.cell_type)
    result["permutations"] = permutations
    result["permutation_seed"] = seed
    result["common_gene_universe_recomputed"] = len(common)
    result["stored_recomputed_max_abs_error"] = maximum_error
    return result


def bubble_size(values: np.ndarray, global_max: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if global_max <= 0:
        return np.full(values.shape, 4.0)
    # Matplotlib scatter uses marker area in points squared.
    return 5.0 + 230.0 * np.clip(values / global_max, 0, 1)


def plot_bubble(
    table: pd.DataFrame,
    output_stem: Path,
    *,
    global_max_similarity: float,
    candidate_only: bool = False,
    correction_label: str = "ComBat v1",
) -> None:
    working = table.copy()
    cell_order = [cell for cell in CELL_ORDER if cell in set(working["cell_type"])]
    protein_order = working["proteomics_module"].drop_duplicates().tolist()
    row_order: list[tuple[str, str]] = []
    boundaries: list[int] = []
    for cell in cell_order:
        cell_rows = working.loc[working["cell_type"].eq(cell)]
        module_order = cell_rows["transcript_module"].drop_duplicates().tolist()
        if candidate_only:
            candidate_modules = set(
                cell_rows.loc[
                    cell_rows["bh_fdr_q_within_figure"].le(0.25),
                    "transcript_module",
                ].astype(str)
            )
            module_order = [module for module in module_order if module in candidate_modules]
        if not module_order:
            continue
        row_order.extend((cell, module) for module in module_order)
        boundaries.append(len(row_order))
    if not row_order:
        return

    row_position = {key: index for index, key in enumerate(row_order)}
    column_position = {module: index for index, module in enumerate(protein_order)}
    x = working["proteomics_module"].map(column_position).to_numpy(dtype=float)
    y = np.asarray(
        [row_position.get((cell, module), np.nan) for cell, module in working[["cell_type", "transcript_module"]].itertuples(index=False, name=None)],
        dtype=float,
    )
    keep = np.isfinite(y)
    x = x[keep]
    y = y[keep]
    values = working.loc[keep, "kme_weighted_jaccard"].to_numpy(dtype=float)
    q_values = working.loc[keep, "bh_fdr_q_within_figure"].to_numpy(dtype=float)

    height = max(8.0, 2.4 + 0.225 * len(row_order))
    figure = plt.figure(figsize=(17.5, height), facecolor="white")
    # Reserve a fixed physical header band.  The lower part of this band is
    # intentionally left for the rotated top-axis module labels.
    header_height = min(0.42, 3.0 / height)
    bottom = min(0.055, 0.42 / height)
    axis_top = 1.0 - header_height
    axis = figure.add_axes([0.31, bottom, 0.675, axis_top - bottom - 0.012])
    header = figure.add_axes([0.025, axis_top, 0.96, header_height])
    header.set_axis_off()
    scatter = axis.scatter(
        x,
        y,
        s=bubble_size(values, global_max_similarity),
        c=q_values,
        cmap="Blues_r",
        norm=Normalize(0, 1),
        edgecolors="#8ea1b2",
        linewidths=0.45,
        alpha=0.98,
    )
    axis.set_xlim(-0.5, len(protein_order) - 0.5)
    axis.set_ylim(len(row_order) - 0.5, -0.5)
    axis.set_xticks(range(len(protein_order)), protein_order, rotation=90, fontsize=8)
    axis.set_yticks(
        range(len(row_order)),
        [f"{CELL_LABELS.get(cell, cell)} | {module}" for cell, module in row_order],
        fontsize=7,
    )
    # Keep the long, rotated protein-module labels below the matrix.  Placing
    # them at the top makes the axis label collide with the size/FDR legends,
    # especially in the compact candidate-only figure.
    axis.xaxis.tick_bottom()
    axis.xaxis.set_label_position("bottom")
    axis.set_xlabel("Proteomics modules", fontsize=10, labelpad=12)
    axis.tick_params(axis="both", length=0)
    axis.grid(False)
    for boundary in boundaries[:-1]:
        axis.axhline(boundary - 0.5, color="#c7d0d8", linewidth=0.75, zorder=0)
    for spine in axis.spines.values():
        spine.set_visible(False)

    threshold_label = "5-cell minimum" if "threshold_5" in output_stem.name else "10-cell minimum"
    title_suffix = "candidate rows (q ≤ 0.25)" if candidate_only else "all module pairs"
    header.text(
        0.0,
        0.93,
        f"CSF transcriptome × protein module pairs — {title_suffix}",
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="#142438",
    )
    header.text(
        0.0,
        0.70,
        f"{correction_label} | {threshold_label} | HVG top 9,000 | Ng_MOG1 and Ng_HD3 removed | other excluded | "
        f"{len(working):,} pair tests across {working['cell_type'].nunique()} cell types\n"
        "Dot area: |kME|-weighted Jaccard   |   Fill: permutation BH-FDR q (darker = smaller q)",
        ha="left",
        va="top",
        fontsize=9.5,
        color="#526174",
    )

    colorbar_axis = header.inset_axes([0.50, 0.31, 0.49, 0.10])
    colorbar = figure.colorbar(
        scatter,
        cax=colorbar_axis,
        orientation="horizontal",
    )
    colorbar.set_label("Permutation BH-FDR q", fontsize=9)
    colorbar.set_ticks([0, 0.05, 0.10, 0.25, 1.0])
    colorbar.ax.tick_params(labelsize=8, length=0)
    colorbar.outline.set_visible(False)

    legend_values = [value for value in (0.05, 0.10, 0.15) if value <= global_max_similarity * 1.05]
    handles = [
        axis.scatter(
            [],
            [],
            s=float(bubble_size(np.asarray([value]), global_max_similarity)[0]),
            facecolor="#718096",
            edgecolor="none",
            label=f"{value:.2f}",
        )
        for value in legend_values
    ]
    if handles:
        header.legend(
            handles=handles,
            title="Weighted Jaccard",
            loc="lower left",
            bbox_to_anchor=(0.18, 0.26),
            ncol=len(handles),
            frameon=False,
            fontsize=8,
            title_fontsize=8,
            handletextpad=0.4,
            columnspacing=1.3,
        )

    header.text(
        0.99,
        0.02,
        "One-sided gene-label permutation; p=(b+1)/(B+1); BH correction within this threshold figure.",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#607080",
    )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in (
        (".png", {"dpi": 300}),
        (".pdf", {}),
        (".svg", {}),
    ):
        figure.savefig(
            output_stem.with_suffix(suffix),
            bbox_inches="tight",
            facecolor="white",
            **kwargs,
        )
    plt.close(figure)


def attach_diagnosis_evidence(
    candidates: pd.DataFrame,
    networks: dict[tuple[str, str], NetworkFiles],
    protein_diagnosis_path: Path,
) -> pd.DataFrame:
    protein = pd.read_csv(protein_diagnosis_path)
    rows: list[dict[str, Any]] = []
    for (threshold, cell_type), frame in candidates.groupby(
        ["threshold", "cell_type"], sort=False
    ):
        network = networks[(str(threshold), str(cell_type))]
        if network.diagnosis_table is None or not network.diagnosis_table.exists():
            continue
        rna = pd.read_csv(network.diagnosis_table)
        for pair in frame.itertuples(index=False):
            rna_module = rna.loc[rna["module_color"].eq(pair.transcript_module)]
            protein_module = protein.loc[
                protein["module_color"].eq(pair.proteomics_module)
            ]
            contrasts = sorted(
                set(rna_module["contrast"].astype(str))
                & set(protein_module["contrast"].astype(str))
            )
            for contrast in contrasts:
                left = rna_module.loc[rna_module["contrast"].eq(contrast)].iloc[0]
                right = protein_module.loc[
                    protein_module["contrast"].eq(contrast)
                ].iloc[0]
                left_r = float(left["point_biserial_r"])
                right_r = float(right["point_biserial_r"])
                rows.append(
                    {
                        "threshold": threshold,
                        "cell_type": cell_type,
                        "transcript_module": pair.transcript_module,
                        "proteomics_module": pair.proteomics_module,
                        "kme_weighted_jaccard": pair.kme_weighted_jaccard,
                        "permutation_p_value": pair.permutation_p_value,
                        "bh_fdr_q_within_figure": pair.bh_fdr_q_within_figure,
                        "contrast": contrast,
                        "rna_point_biserial_r": left_r,
                        "rna_diagnosis_q": float(left["bh_fdr_q_within_network"]),
                        "protein_point_biserial_r": right_r,
                        "protein_diagnosis_q": float(
                            right["bh_fdr_q_within_network"]
                        ),
                        "same_direction": bool(np.sign(left_r) == np.sign(right_r)),
                        "both_diagnosis_q_le_0p10": bool(
                            float(left["bh_fdr_q_within_network"]) <= 0.10
                            and float(right["bh_fdr_q_within_network"]) <= 0.10
                        ),
                    }
                )
    return pd.DataFrame(rows)


def write_candidate_shared_genes(
    candidates: pd.DataFrame,
    networks: dict[tuple[str, str], NetworkFiles],
    output_path: Path,
) -> None:
    collected: list[pd.DataFrame] = []
    for (threshold, cell_type), frame in candidates.groupby(
        ["threshold", "cell_type"], sort=False
    ):
        network = networks[(str(threshold), str(cell_type))]
        contributions = pd.read_csv(network.contribution_table)
        selected = contributions.merge(
            frame[
                [
                    "transcript_module",
                    "proteomics_module",
                    "kme_weighted_jaccard",
                    "permutation_p_value",
                    "bh_fdr_q_within_figure",
                ]
            ],
            on=["transcript_module", "proteomics_module"],
            how="inner",
            validate="many_to_one",
        )
        if selected.empty:
            continue
        selected.insert(0, "threshold", threshold)
        selected.insert(1, "cell_type", cell_type)
        selected["shared_gene_rank_within_pair"] = selected.groupby(
            ["transcript_module", "proteomics_module"]
        )["weighted_intersection"].rank(method="first", ascending=False)
        collected.append(selected)
    result = pd.concat(collected, ignore_index=True) if collected else pd.DataFrame()
    result.to_csv(
        output_path,
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )


def bubble_stem(output_dir: Path, threshold: str) -> Path:
    return output_dir / (
        f"CSF_{threshold}_hvg_top9000_after_outlier_"
        "kME_weighted_Jaccard_permutation_bubble"
    )


def write_viewer_statistics(output_dir: Path, batch_correction: str) -> None:
    """Package saved p/q values for file:// use, keyed by overlap record and pair."""
    thresholds = {}
    for threshold in THRESHOLDS:
        stem = bubble_stem(output_dir, threshold)
        stats_path = stem.with_name(stem.name + "_statistics.csv.gz")
        table = pd.read_csv(stats_path)
        columns = [
            "record_id", "transcript_module", "proteomics_module",
            "permutation_p_value", "bh_fdr_q_within_figure",
        ]
        keys = columns[:3]
        if table.empty or table.duplicated(keys).any():
            raise ValueError(f"Empty or duplicate module pairs in {stats_path}")
        is_v2 = table["record_id"].str.contains("::combat_v2_cell_diagnosis_tissue::", regex=False)
        if not is_v2.eq(batch_correction == "combat_v2_cell_diagnosis_tissue").all():
            raise ValueError(f"Unexpected batch correction in {stats_path}")
        thresholds[threshold] = {
            "module_pairs": len(table),
            "q_le_0p25": int(table["bh_fdr_q_within_figure"].le(0.25).sum()),
            "stem": stem.name,
            "files": sorted(path.name for path in output_dir.glob(stem.name + "*")),
            "pairs": [
                [record_id + "::after_outlier", rna, protein, p, q]
                for record_id, rna, protein, p, q in table[columns].itertuples(index=False, name=None)
            ],
        }
    (output_dir / "viewer-statistics.js").write_text(
        "window.FINAL_JACCARD_STATISTICS=window.FINAL_JACCARD_STATISTICS||{};\n"
        f"window.FINAL_JACCARD_STATISTICS[{json.dumps(batch_correction)}]="
        + json.dumps(thresholds, separators=(",", ":"), allow_nan=False)
        + ";\n",
        encoding="utf-8",
    )


def run(
    *, permutations: int, seed: int, output_dir: Path,
    batch_correction: str = DEFAULT_CORRECTION, dashboard_root: Path = ROOT,
) -> None:
    records = read_dashboard_records(dashboard_root)
    networks = select_networks(records, batch_correction, dashboard_root)
    if {network.threshold for network in networks} != set(THRESHOLDS):
        raise RuntimeError(f"Missing after-outlier CSF HVG-9000 networks for {batch_correction}")
    proteomics_record = next(
        record for record in records if record["id"] == PROTEOMICS_RECORD_ID
    )
    protein_files = file_map(proteomics_record, dashboard_root)
    protein_source = pd.read_csv(protein_files["official_kme"])
    protein_collapsed, protein_universe, collapse_audit = collapse_proteomics_to_gene(
        protein_source
    )
    first_pairs = pd.read_csv(networks[0].pair_table)
    protein_module_order = (
        first_pairs["proteomics_module"].drop_duplicates().astype(str).tolist()
    )

    results: list[pd.DataFrame] = []
    for index, network in enumerate(networks):
        print(
            f"[{index + 1:02d}/{len(networks):02d}] {network.threshold} / "
            f"{network.cell_type}: {permutations:,} permutations",
            flush=True,
        )
        results.append(
            analyse_network(
                network,
                protein_collapsed,
                protein_universe,
                protein_module_order,
                permutations=permutations,
                seed=seed + index,
            )
        )
    combined = pd.concat(results, ignore_index=True)
    combined["bh_fdr_q_within_figure"] = np.nan
    for threshold, index in combined.groupby("threshold", sort=False).groups.items():
        combined.loc[index, "bh_fdr_q_within_figure"] = benjamini_hochberg(
            combined.loc[index, "permutation_p_value"].to_numpy(dtype=float)
        )
    combined["similarity_q_le_0p05"] = combined[
        "bh_fdr_q_within_figure"
    ].le(0.05)
    combined["similarity_q_le_0p10"] = combined[
        "bh_fdr_q_within_figure"
    ].le(0.10)
    combined["similarity_q_le_0p25"] = combined[
        "bh_fdr_q_within_figure"
    ].le(0.25)
    global_max = float(combined["kme_weighted_jaccard"].max())
    output_dir.mkdir(parents=True, exist_ok=True)

    network_lookup = {
        (network.threshold, network.cell_type): network for network in networks
    }
    summary: dict[str, Any] = {
        "analysis": "CSF cross-omics module overlap permutation test",
        "batch_correction": batch_correction,
        "strategy": STRATEGY,
        "sample_exclusion": ["Ng_MOG1", "Ng_HD3"],
        "cell_type_exclusion": ["other"],
        "permutations": permutations,
        "base_seed": seed,
        "p_value": "(number of permuted weighted Jaccard values >= observed + 1) / (B + 1)",
        "null_model": (
            "permute RNA gene identities within each network-specific common-gene universe; "
            "preserve both omics module sizes and absolute own-module kME weights"
        ),
        "multiple_testing": "Benjamini-Hochberg within all cell-type x RNA-module x protein-module pairs in each threshold figure",
        "candidate_flags": [0.05, 0.10, 0.25],
        "protein_module_order": protein_module_order,
        "proteomics_collapse_audit": collapse_audit,
        "references": {
            "permutation_p_value": "Phipson and Smyth 2010, doi:10.2202/1544-6115.1585",
            "fdr": "Benjamini and Hochberg 1995, doi:10.1111/j.2517-6161.1995.tb02031.x",
        },
        "thresholds": {},
    }

    for threshold in THRESHOLDS:
        table = combined.loc[combined["threshold"].eq(threshold)].copy()
        stem = bubble_stem(output_dir, threshold)
        table.to_csv(
            stem.with_name(stem.name + "_statistics.csv.gz"),
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        candidates = table.loc[table["similarity_q_le_0p25"]].sort_values(
            ["bh_fdr_q_within_figure", "kme_weighted_jaccard"],
            ascending=[True, False],
            kind="stable",
        )
        candidates.to_csv(stem.with_name(stem.name + "_candidates_q0p25.csv"), index=False)
        plot_bubble(
            table,
            stem,
            global_max_similarity=global_max,
            candidate_only=False,
            correction_label=CORRECTION_LABELS[batch_correction],
        )
        plot_bubble(
            table,
            stem.with_name(stem.name + "_candidate_rows_q0p25"),
            global_max_similarity=global_max,
            candidate_only=True,
            correction_label=CORRECTION_LABELS[batch_correction],
        )
        write_candidate_shared_genes(
            candidates,
            network_lookup,
            stem.with_name(stem.name + "_candidate_shared_genes.csv.gz"),
        )
        diagnosis = attach_diagnosis_evidence(
            candidates,
            network_lookup,
            protein_files["diagnosis_pairwise_statistics"],
        )
        diagnosis.to_csv(
            stem.with_name(stem.name + "_candidate_diagnosis_concordance.csv"),
            index=False,
        )
        summary["thresholds"][threshold] = {
            "cell_types": int(table["cell_type"].nunique()),
            "module_pairs": int(len(table)),
            "maximum_weighted_jaccard": float(
                table["kme_weighted_jaccard"].max()
            ),
            "minimum_permutation_p": float(table["permutation_p_value"].min()),
            "minimum_bh_fdr_q": float(table["bh_fdr_q_within_figure"].min()),
            "q_le_0p05": int(table["similarity_q_le_0p05"].sum()),
            "q_le_0p10": int(table["similarity_q_le_0p10"].sum()),
            "q_le_0p25": int(table["similarity_q_le_0p25"].sum()),
        }
    (output_dir / "analysis_method_and_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_viewer_statistics(output_dir, batch_correction)
    print(json.dumps(summary["thresholds"], indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--dashboard-root", type=Path, default=ROOT)
    parser.add_argument("--batch-correction", choices=CORRECTION_LABELS, default=DEFAULT_CORRECTION)
    parser.add_argument("--from-saved-statistics", action="store_true")
    args = parser.parse_args()
    if args.from_saved_statistics:
        write_viewer_statistics(args.output_dir.resolve(), args.batch_correction)
    else:
        run(
            permutations=args.permutations,
            seed=args.seed,
            output_dir=args.output_dir.resolve(),
            batch_correction=args.batch_correction,
            dashboard_root=args.dashboard_root.resolve(),
        )
