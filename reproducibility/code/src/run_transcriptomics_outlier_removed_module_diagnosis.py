"""Run kME and diagnosis downstream for after-outlier findModules results."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
from typing import Any

import pandas as pd

from run_official_kme_module_diagnosis import (
    CONDITION,
    run_network as run_official_network,
)
from run_pairwise_hd_module_diagnosis import run_network as run_pairwise_network


def network_root(stage: Path, row: Any) -> Path:
    root = stage / row.threshold / row.strategy / row.mode / row.specimen
    if row.mode == "cell_separate":
        root /= row.cell_type
    return root


def registry(stage: Path) -> pd.DataFrame:
    frame = pd.read_csv(stage / "outlier_removal_find_modules_summary.csv")
    strategy_rank = {
        "03_hvg_top5000": 0,
        "02_hvg_top9000": 1,
        "01_cpm_detection": 2,
    }
    frame = frame.loc[frame["status"].eq("complete_findModules")].copy()
    frame["_strategy"] = frame["strategy"].map(strategy_rank)
    frame["_threshold"] = frame["threshold"].map(
        {"threshold_5_cells": 0, "threshold_10_cells": 1}
    )
    return frame.sort_values(
        ["_strategy", "_threshold", "mode", "specimen", "cell_type"],
        kind="stable",
    ).reset_index(drop=True)


def dimensions(row: Any) -> dict[str, str]:
    return {
        "threshold": str(row.threshold),
        "strategy": str(row.strategy),
        "mode": str(row.mode),
        "specimen": str(row.specimen),
        "cell_type": str(row.cell_type),
        "condition": CONDITION,
        "sample_exclusion": "Ng_MOG1;Ng_HD3",
        "sample_tree_cut": "none_cut_inf",
    }


def process(run_root: Path, stage: Path, row: Any, overwrite: bool) -> dict[str, Any]:
    root = network_root(stage, row)
    branch = root / "outlier_removal"
    module_root = branch / "05_find_modules" / CONDITION
    output_root = branch / "06_kme_module_diagnosis" / CONDITION
    matrix_path = branch / "03_preprocess/01_cut_diagnostic/matrix_before_cut.csv.gz"
    sample_info_path = module_root / "sample_info.csv"
    values = dimensions(row)
    official = run_official_network(
        run_root=run_root,
        omics="transcriptomics",
        network_id=str(row.network_id),
        source_model=None,
        source_matrix=matrix_path,
        sample_info_path=sample_info_path,
        assignments_path=module_root / "module_assignments.csv.gz",
        assignment_feature_column="gene_id",
        eigengenes_path=module_root / "module_eigengenes.csv",
        output_root=output_root,
        dimensions=values,
        overwrite=overwrite,
    )
    pairwise = run_pairwise_network(
        run_root=run_root,
        omics="transcriptomics",
        network_id=str(row.network_id),
        source_model=None,
        source_matrix=matrix_path,
        sample_info_path=sample_info_path,
        eigengenes_path=module_root / "module_eigengenes.csv",
        output_root=output_root,
        dimensions=values,
        overwrite=overwrite,
    )
    return {
        "threshold": row.threshold,
        "strategy": row.strategy,
        "mode": row.mode,
        "specimen": row.specimen,
        "cell_type": row.cell_type,
        "network_id": row.network_id,
        "official_status": official.get("status"),
        "pairwise_status": pairwise.get("status"),
        "samples": official.get("samples"),
        "features": official.get("features"),
        "module_eigengenes": official.get("module_eigengenes"),
    }


def collect(stage: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in registry(stage).itertuples(index=False):
        root = network_root(stage, row) / "outlier_removal/06_kme_module_diagnosis" / CONDITION
        official_path = root / "analysis_summary.json"
        pairwise_path = root / "pairwise_hd_analysis_summary.json"
        official = json.loads(official_path.read_text(encoding="utf-8")) if official_path.exists() else {}
        pairwise = json.loads(pairwise_path.read_text(encoding="utf-8")) if pairwise_path.exists() else {}
        rows.append(
            {
                "threshold": row.threshold,
                "strategy": row.strategy,
                "mode": row.mode,
                "specimen": row.specimen,
                "cell_type": row.cell_type,
                "network_id": row.network_id,
                "official_status": official.get("status", "missing"),
                "pairwise_status": pairwise.get("status", "missing"),
                "samples": official.get("samples"),
                "features": official.get("features"),
                "module_eigengenes": official.get("module_eigengenes"),
            }
        )
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> None:
    run_root = args.run_root.resolve()
    stage = run_root / "experiments/04_variance_filtered/transcriptomics"
    if args.summary_only:
        frame = collect(stage)
        path = stage / "outlier_removal_module_diagnosis_summary.csv"
        frame.to_csv(path, index=False)
        print(f"Updated {path} ({len(frame)} rows)", flush=True)
        return
    frame = registry(stage)
    if args.only_network:
        frame = frame.loc[frame["network_id"].eq(args.only_network)].copy()
    if args.limit is not None:
        frame = frame.head(args.limit).copy()
    if args.worker_count < 1 or not 0 <= args.worker_index < args.worker_count:
        raise ValueError("worker index/count are invalid")
    frame = frame.iloc[args.worker_index :: args.worker_count]
    for position, row in enumerate(frame.itertuples(index=False), start=1):
        print(
            f"[after-outlier downstream {position}/{len(frame)}] {row.strategy} / {row.network_id}",
            flush=True,
        )
        try:
            result = process(run_root, stage, row, args.overwrite)
            print(
                f"  official={result['official_status']} | pairwise={result['pairwise_status']}",
                flush=True,
            )
        except BaseException as error:
            root = network_root(stage, row) / "outlier_removal/06_kme_module_diagnosis" / CONDITION
            root.mkdir(parents=True, exist_ok=True)
            (root / "downstream_error.json").write_text(
                json.dumps(
                    {"status": "failed", "error_type": type(error).__name__, "error": str(error)},
                    indent=2,
                ),
                encoding="utf-8",
            )
            print(f"  FAILED: {type(error).__name__}: {error}", flush=True)
        gc.collect()
    print(
        f"Completed worker {args.worker_index + 1}/{args.worker_count}; run --summary-only after all workers.",
        flush=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    parser.add_argument("--only-network")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
