"""Add CSF-RNA versus PB-protein comparisons to the offline Final Jaccard view.

The established cross-omics pipeline compares specimen-matched networks.  The
Final Jaccard page, however, needs a controlled sensitivity view in which the
same CSF transcriptomics modules are compared with either the CSF or the PB
proteomics network.  This script computes only the missing CSF-RNA x PB-protein
branch with the exact same absolute-own-module-kME weighted Jaccard definition,
writes auditable CSV results, and registers file:// compatible dashboard shards.

Existing CSF-RNA x CSF-protein records are never recalculated or overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from run_cross_omics_kme_weighted_jaccard import (
    PROTEOMICS_REFERENCE_CONDITION,
    TRANSCRIPT_CONDITION,
    collapse_proteomics_to_gene,
    compare_modules,
    module_weights,
    normalize_gene,
)


FINAL_THRESHOLDS = {"threshold_5_cells", "threshold_10_cells"}
FINAL_STRATEGIES = {"01_cpm_detection", "02_hvg_top9000", "03_hvg_top5000"}
TRANSCRIPT_SPECIMEN = "CSF"
PROTEIN_SPECIMEN = "PB"
PAIR_COLUMNS = [
    "transcript_module",
    "proteomics_module",
    "transcript_module_genes",
    "proteomics_module_genes",
    "shared_genes",
    "jaccard",
    "weighted_intersection",
    "weighted_union",
    "weighted_jaccard",
    "transcript_best",
    "proteomics_best",
    "reciprocal_best",
]


def _read_js_assignment(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        payload = text.split("=", 1)[1].strip()
    except IndexError as exc:
        raise ValueError(f"Not a JavaScript assignment: {path}") from exc
    if payload.endswith(";"):
        payload = payload[:-1]
    return json.loads(payload)


def _finite(value: Any) -> float | None:
    if value is None or pd.isna(value):
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def _base_record_id(row: Any) -> str:
    return "::".join(
        [
            "transcriptomics_variance",
            str(row.threshold),
            str(row.strategy),
            str(row.mode),
            str(row.specimen),
            str(row.cell_type),
        ]
    )


def _network_root(stage_root: Path, row: Any, after_outlier: bool) -> Path:
    root = stage_root / row.threshold / row.strategy / row.mode / row.specimen
    if row.mode == "cell_separate":
        root = root / row.cell_type
    return root / "outlier_removal" if after_outlier else root


def _result_root(run_root: Path, row: Any, branch: str) -> Path:
    return (
        run_root
        / "experiments/05_cross_omics_module_similarity/final_jaccard_protein_scope"
        / branch
        / f"transcript_{TRANSCRIPT_SPECIMEN}"
        / f"protein_{PROTEIN_SPECIMEN}"
        / row.threshold
        / row.strategy
        / row.cell_type
    )


def _pair_rows(pairs: pd.DataFrame) -> list[list[Any]]:
    ordered = pairs.sort_values(
        ["kme_weighted_jaccard", "transcript_module", "proteomics_module"],
        ascending=[False, True, True],
        kind="stable",
    )
    rows: list[list[Any]] = []
    for row in ordered.itertuples(index=False):
        rows.append(
            [
                str(row.transcript_module),
                str(row.proteomics_module),
                int(row.transcript_module_genes_in_common),
                int(row.proteomics_module_genes_in_common),
                int(row.shared_genes),
                _finite(row.jaccard),
                _finite(row.weighted_intersection),
                _finite(row.weighted_union),
                _finite(row.kme_weighted_jaccard),
                bool(row.transcript_best),
                bool(row.proteomics_best),
                bool(row.reciprocal_best),
            ]
        )
    return rows


def _contribution_rows(table: pd.DataFrame) -> list[list[Any]]:
    if table.empty:
        return []
    ordered = table.sort_values(
        ["transcript_module", "proteomics_module", "weighted_intersection", "gene"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    return [
        [
            str(row.transcript_module),
            str(row.proteomics_module),
            str(row.gene),
            _finite(row.transcript_abs_kME),
            _finite(row.proteomics_abs_kME),
            _finite(row.weighted_intersection),
            _finite(row.weighted_union),
        ]
        for row in ordered.itertuples(index=False)
    ]


def _annotate_existing_records(index: dict[str, Any]) -> None:
    """Make the protein network explicit for legacy specimen-matched records."""

    for record_id, entry in index.get("records", {}).items():
        entry.setdefault("protein_specimen", entry.get("specimen"))
        entry.setdefault("transcript_record_id", record_id)
        parts = str(record_id).split("::")
        after = "after_outlier" in parts
        base = [part for part in parts if part != "after_outlier"]
        if len(base) >= 6 and base[0] == "transcriptomics_variance":
            entry.setdefault("threshold", base[1])
            entry.setdefault("strategy", base[2])
            entry.setdefault("mode", base[3])
            entry.setdefault("cell_type", base[5])
        entry.setdefault(
            "analysis_branch", "after_outlier_removal" if after else "baseline"
        )


def build(run_root: Path, dashboard_root: Path) -> dict[str, Any]:
    transcript_stage = run_root / "experiments/04_variance_filtered/transcriptomics"
    protein_path = (
        run_root
        / "experiments/01_batch_corrected/proteomics/06_kme_module_diagnosis"
        / PROTEIN_SPECIMEN
        / PROTEOMICS_REFERENCE_CONDITION
        / "feature_module_kme.csv.gz"
    )
    overlap_root = dashboard_root / "cross-omics-overlap"
    shard_root = overlap_root / "shards"
    index_path = overlap_root / "overlap-index.js"
    if not index_path.is_file():
        raise FileNotFoundError(index_path)
    if not protein_path.is_file():
        raise FileNotFoundError(protein_path)

    protein_table = pd.read_csv(protein_path)
    protein_collapsed, protein_universe, collapse_audit = collapse_proteomics_to_gene(
        protein_table
    )
    protein_collapsed["is_raw_unassigned"] = False
    protein_module_names = list(
        dict.fromkeys(
            protein_table.loc[
                ~protein_table["is_raw_unassigned"].astype(bool)
                & protein_table["official_assigned_abs_kME"].notna(),
                "module_color",
            ]
            .astype(str)
            .tolist()
        )
    )

    index = _read_js_assignment(index_path)
    _annotate_existing_records(index)
    shard_root.mkdir(parents=True, exist_ok=True)
    generated_ids: set[str] = set()
    summary_rows: list[dict[str, Any]] = []

    for branch, registry_name in (
        ("baseline", "find_modules_summary.csv"),
        ("after_outlier_removal", "outlier_removal_find_modules_summary.csv"),
    ):
        after_outlier = branch == "after_outlier_removal"
        registry = pd.read_csv(transcript_stage / registry_name)
        registry = registry.loc[
            registry["status"].eq("complete_findModules")
            & registry["specimen"].eq(TRANSCRIPT_SPECIMEN)
            & registry["mode"].eq("cell_separate")
            & registry["threshold"].isin(FINAL_THRESHOLDS)
            & registry["strategy"].isin(FINAL_STRATEGIES)
            & ~registry["cell_type"].astype(str).str.lower().eq("other")
        ].copy()

        for row in registry.itertuples(index=False):
            transcript_path = (
                _network_root(transcript_stage, row, after_outlier)
                / "06_kme_module_diagnosis"
                / TRANSCRIPT_CONDITION
                / "feature_module_kme.csv.gz"
            )
            if not transcript_path.is_file():
                continue
            transcript = pd.read_csv(transcript_path)
            transcript["gene_symbol"] = transcript["gene_id"].map(normalize_gene)
            transcript_universe = set(transcript["gene_symbol"].dropna().astype(str))
            common = transcript_universe & protein_universe
            transcript_module_names = list(
                dict.fromkeys(
                    transcript.loc[
                        ~transcript["is_raw_unassigned"].astype(bool)
                        & transcript["official_assigned_abs_kME"].notna(),
                        "module_color",
                    ]
                    .astype(str)
                    .tolist()
                )
            )
            transcript_modules, transcript_counts = module_weights(
                transcript,
                "gene_symbol",
                "official_assigned_abs_kME",
                common,
                transcript_module_names,
            )
            protein_modules, protein_counts = module_weights(
                protein_collapsed,
                "gene_symbol",
                "official_assigned_abs_kME",
                common,
                protein_module_names,
            )
            pairs, contributions = compare_modules(
                transcript_modules,
                protein_modules,
                transcript_counts,
                protein_counts,
                len(common),
            )
            if pairs.empty:
                continue

            result_root = _result_root(run_root, row, branch)
            result_root.mkdir(parents=True, exist_ok=True)
            pair_path = result_root / "module_pair_similarity.csv"
            contribution_path = result_root / "all_pair_gene_contributions.csv.gz"
            pairs.to_csv(pair_path, index=False)
            contributions.to_csv(
                contribution_path,
                index=False,
                compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
            )

            base_record_id = _base_record_id(row)
            transcript_record_id = (
                f"{base_record_id}::after_outlier" if after_outlier else base_record_id
            )
            record_id = f"{transcript_record_id}::protein_{PROTEIN_SPECIMEN}"
            filename = hashlib.sha256(record_id.encode("utf-8")).hexdigest()[:24] + ".js"
            shard_path = shard_root / filename
            payload = {
                "columns": [
                    "transcript_module",
                    "proteomics_module",
                    "gene",
                    "transcript_abs_kme",
                    "proteomics_abs_kme",
                    "weighted_intersection",
                    "weighted_union",
                ],
                "rows": _contribution_rows(contributions),
            }
            shard_path.write_text(
                "window.CROSS_OMICS_OVERLAP_SHARDS=window.CROSS_OMICS_OVERLAP_SHARDS||{};"
                f"window.CROSS_OMICS_OVERLAP_SHARDS[{json.dumps(record_id)}]="
                + json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
                + ";\n",
                encoding="utf-8",
            )

            index["records"][record_id] = {
                "title": (
                    f"{row.cell_type} - CSF RNA x PB protein"
                    + (" - After outlier removal" if after_outlier else "")
                ),
                "network_id": row.network_id,
                "specimen": TRANSCRIPT_SPECIMEN,
                "protein_specimen": PROTEIN_SPECIMEN,
                "transcript_record_id": transcript_record_id,
                "cell_type": row.cell_type,
                "threshold": row.threshold,
                "strategy": row.strategy,
                "mode": row.mode,
                "analysis_branch": branch,
                "proteomics_condition": PROTEOMICS_REFERENCE_CONDITION,
                "path": f"cross-omics-overlap/shards/{filename}",
                "transcript_modules": transcript_module_names,
                "proteomics_modules": protein_module_names,
                "pair_columns": PAIR_COLUMNS,
                "pairs": _pair_rows(pairs),
                "overlap_rows": int(len(contributions)),
            }
            generated_ids.add(record_id)

            summary = {
                "record_id": record_id,
                "transcript_record_id": transcript_record_id,
                "analysis_branch": branch,
                "threshold": row.threshold,
                "strategy": row.strategy,
                "cell_type": row.cell_type,
                "transcript_specimen": TRANSCRIPT_SPECIMEN,
                "protein_specimen": PROTEIN_SPECIMEN,
                "common_gene_universe": len(common),
                "transcript_module_count": len(transcript_module_names),
                "protein_module_count": len(protein_module_names),
                "module_pairs": len(pairs),
                "overlap_rows": len(contributions),
                "maximum_weighted_jaccard": _finite(
                    pairs["kme_weighted_jaccard"].max()
                ),
                "pair_table": pair_path.as_posix(),
                "contribution_table": contribution_path.as_posix(),
                "transcript_kme": transcript_path.as_posix(),
                "protein_kme": protein_path.as_posix(),
            }
            summary_rows.append(summary)
            (result_root / "analysis_summary.json").write_text(
                json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )

    # Remove only supplemental records that this builder owns but no longer emits.
    stale_ids = [
        record_id
        for record_id in index.get("records", {})
        if record_id.endswith(f"::protein_{PROTEIN_SPECIMEN}")
        and record_id not in generated_ids
    ]
    for record_id in stale_ids:
        entry = index["records"].pop(record_id)
        stale_path = dashboard_root / str(entry.get("path", ""))
        if stale_path.is_file():
            stale_path.unlink()

    index["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    index_path.write_text(
        "window.CROSS_OMICS_OVERLAP_INDEX="
        + json.dumps(index, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )

    dashboard_data_path = dashboard_root / "dashboard-data.js"
    if dashboard_data_path.is_file():
        dashboard_data = _read_js_assignment(dashboard_data_path)
        dashboard_data["cross_omics_overlap_explorer"] = {
            "records": len(index["records"]),
            "pairs": sum(len(entry.get("pairs", [])) for entry in index["records"].values()),
            "overlap_rows": sum(
                int(entry.get("overlap_rows", 0) or 0)
                for entry in index["records"].values()
            ),
            "index": "cross-omics-overlap/overlap-index.js",
        }
        dashboard_data_path.write_text(
            "window.EXPERIMENT_DASHBOARD_DATA = "
            + json.dumps(dashboard_data, ensure_ascii=False, separators=(",", ":"))
            + ";\n",
            encoding="utf-8",
        )

    output_root = (
        run_root
        / "experiments/05_cross_omics_module_similarity/final_jaccard_protein_scope"
    )
    summary_table = pd.DataFrame(summary_rows)
    summary_path = output_root / "CSF_RNA_PB_protein_comparison_summary.csv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_table.to_csv(summary_path, index=False)
    method = {
        "method": "absolute-own-module-kME weighted Jaccard",
        "formula": "sum_g min(w_RNA,g, w_protein,g) / sum_g max(w_RNA,g, w_protein,g)",
        "transcript_specimen": TRANSCRIPT_SPECIMEN,
        "protein_specimen": PROTEIN_SPECIMEN,
        "proteomics_reference_condition": PROTEOMICS_REFERENCE_CONDITION,
        "universe": "per-network intersection of measurable RNA genes and simple protein gene symbols",
        "protein_probe_collapse": collapse_audit.get("selection_rule"),
        "purpose": "protein-specimen sensitivity comparison for the Final Jaccard dashboard",
        "records": len(summary_rows),
    }
    (output_root / "method.json").write_text(
        json.dumps(method, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "records": len(summary_rows),
        "baseline_records": sum(
            row["analysis_branch"] == "baseline" for row in summary_rows
        ),
        "after_outlier_records": sum(
            row["analysis_branch"] == "after_outlier_removal"
            for row in summary_rows
        ),
        "index": index_path.as_posix(),
        "summary": summary_path.as_posix(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--dashboard-root", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build(args.run_root.resolve(), args.dashboard_root.resolve())
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
