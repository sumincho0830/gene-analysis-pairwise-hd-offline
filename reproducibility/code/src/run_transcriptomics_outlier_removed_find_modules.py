"""Run findModules from the cut=inf preprocess matrix after two exclusions.

The input is the audited ``matrix_before_cut.csv.gz`` produced after removing
Ng_MOG1 and Ng_HD3 and running PyWGCNA preprocess with ``cut=inf``.  No finite
sample-tree cut, extra feature filtering, value transform, scaling, or
imputation is introduced here.  Results are additive under each network's
``outlier_removal`` branch so the earlier variance-filtered results remain
unchanged.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import html
import json
from pathlib import Path
import time
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_provenance import write_run_manifest
from run_transcriptomics_stage1 import PipelineConfig
from run_transcriptomics_variance_filtered import (
    STRATEGIES,
    Network,
    load_networks,
    network_metadata,
)
from run_transcriptomics_variance_find_modules import (
    CONDITION,
    export_tables,
    generate_dynamic_me_dendrogram,
    plot_gene_dendrogram_modules,
    plot_module_sizes,
    plot_soft_threshold,
    run_find_modules,
    validate_input_model,
)
from transcriptomics_qc import build_source_cohort_map, read_expression_matrix


STRATEGY_ORDER = (
    "03_hvg_top5000",
    "02_hvg_top9000",
    "01_cpm_detection",
)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def strategy_for(network: Network) -> str:
    matches = [value for value in STRATEGIES if value in network.output_root.parts]
    if len(matches) != 1:
        raise ValueError(f"Cannot resolve strategy from {network.output_root}")
    return matches[0]


def ordered_networks(config: PipelineConfig) -> list[tuple[Network, str]]:
    networks = [(network, strategy_for(network)) for network in load_networks(config, STRATEGIES)]
    strategy_rank = {value: index for index, value in enumerate(STRATEGY_ORDER)}
    return sorted(
        networks,
        key=lambda item: (
            strategy_rank[item[1]],
            0 if item[0].threshold == "threshold_5_cells" else 1,
            0 if item[0].mode == "cell_separate" else 1,
            item[0].specimen,
            item[0].cell_type.lower(),
        ),
    )


def compact(network: Network, strategy: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "threshold": network.threshold,
        "strategy": strategy,
        "mode": network.mode,
        "specimen": network.specimen,
        "cell_type": network.cell_type,
        "network_id": network.network_id,
        **payload,
    }


def output_root(network: Network) -> Path:
    return network.output_root / "outlier_removal/05_find_modules" / CONDITION


def process_network(
    run_root: Path,
    config: PipelineConfig,
    network: Network,
    strategy: str,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
) -> dict[str, Any]:
    root = output_root(network)
    summary_path = root / "analysis_summary.json"
    manifest_path = root / "run_manifest.json"
    if summary_path.exists() and manifest_path.exists():
        previous = json.loads(summary_path.read_text(encoding="utf-8"))
        # Only a completed module result is immutable.  A previous skip can
        # become stale when the upstream preprocess status is finalized, so it
        # must be re-evaluated against preprocess_summary.json on every run.
        if previous.get("status") == "complete_findModules":
            return compact(network, strategy, previous)

    diagnostic = network.output_root / "outlier_removal/03_preprocess/01_cut_diagnostic"
    preprocess_path = diagnostic / "preprocess_summary.json"
    preprocess = json.loads(preprocess_path.read_text(encoding="utf-8"))
    root.mkdir(parents=True, exist_ok=True)
    (root / "figures").mkdir(parents=True, exist_ok=True)
    if preprocess.get("status") not in {
        "complete_preprocess_cut_inf",
        "complete_cut_diagnostic_pending_manual_review",
    }:
        payload = {
            "status": "not_run_preprocess_incomplete",
            "reason": preprocess.get("reason", preprocess.get("status")),
            "condition": CONDITION,
            "findModules_run": False,
            "finite_sample_cut_applied": False,
        }
        write_json(summary_path, payload)
        write_run_manifest(
            manifest_path,
            pipeline_name="transcriptomics_after_outlier_removal_find_modules_skip",
            run_root=run_root,
            parameters={"condition": CONDITION},
            generation_rules={"reason": payload["reason"]},
            input_paths=[preprocess_path],
            output_paths=[summary_path],
            code_paths=[Path(__file__)],
            packages=("PyWGCNA",),
        )
        return compact(network, strategy, payload)

    from PyWGCNA import WGCNA

    matrix_path = diagnostic / "matrix_before_cut.csv.gz"
    matrix = read_expression_matrix(matrix_path)
    metadata = network_metadata(network, matrix, clinical, cohort_map)
    sample_columns = [
        "project_id",
        "patient_id",
        "cohort",
        "specimen",
        "cell_type",
        "diagnosis",
        "sex",
        "age",
    ]
    sample_info = metadata.loc[matrix.index, sample_columns].copy()
    sample_info.index = matrix.index.astype(str)
    sample_info.index.name = "sample_id"
    sample_info_path = root / "sample_info.csv"
    sample_info.to_csv(sample_info_path)

    started = time.perf_counter()
    model = WGCNA(
        name=f"TRANSCRIPTOMICS_{network.network_id}_{strategy}_after_outlier_removal",
        geneExp=matrix,
        sampleInfo=sample_info,
        species="Human",
        save=False,
        outputPath=root.resolve().as_posix() + "/",
        figureType="png",
        cut=float("inf"),
    )
    validate_input_model(model, expected_samples=matrix.shape[0], expected_genes=matrix.shape[1])
    raw_dynamic, fallback_reason = run_find_modules(model)
    assignments, module_summary = export_tables(model, raw_dynamic, root)
    generate_dynamic_me_dendrogram(
        model.datExpr.to_df(),
        assignments,
        root,
        me_diss_threshold=float(model.MEDissThres),
    )

    label = (
        f"after Ng_MOG1/Ng_HD3 removal; {network.threshold} / {strategy} / "
        f"{network.mode} / {network.specimen} / {network.cell_type}"
    )
    figures = root / "figures"
    plot_gene_dendrogram_modules(
        model,
        figures / "gene_dendrogram_modules.png",
        f"Transcriptomics gene dendrogram and modules: {label}",
    )
    plot_soft_threshold(
        model,
        figures / "soft_threshold.png",
        f"Transcriptomics soft-threshold: {label}",
    )
    plot_module_sizes(
        module_summary,
        figures / "module_sizes.png",
        f"Transcriptomics module sizes: {label}",
    )

    proper = module_summary.loc[~module_summary["is_fully_raw_unassigned"]]
    total_assigned = int(proper["assigned_feature_count"].sum())
    median_kme = assignments.loc[
        ~assignments["is_raw_unassigned"], "assigned_abs_kME"
    ].median()
    payload = {
        "status": "complete_findModules",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "condition": CONDITION,
        "samples": int(model.datExpr.shape[0]),
        "genes": int(model.datExpr.shape[1]),
        "confirmed_subjects_removed": ["Ng_MOG1", "Ng_HD3"],
        "finite_sample_cut_applied": False,
        "preprocess_input": "cut=inf matrix_before_cut",
        "feature_filtering_before_findModules": "none",
        "expression_transform_before_findModules": "none",
        "selected_power": int(model.power),
        "tested_powers": [int(value) for value in model.powers],
        "RsquaredCut": float(model.RsquaredCut),
        "MeanCut": float(model.MeanCut),
        "networkType": str(model.networkType),
        "TOMType": str(model.TOMType),
        "deepSplit": 2,
        "minModuleSize": int(model.minModuleSize),
        "MEDissThres": float(model.MEDissThres),
        "pamRespectsDendro": False,
        "dynamic_module_count": int(
            assignments.loc[~assignments["is_raw_unassigned"], "raw_dynamic_label"].nunique()
        ),
        "final_module_count": int(len(proper)),
        "raw_unassigned_gene_count": int(assignments["is_raw_unassigned"].sum()),
        "raw_unassigned_fraction": float(assignments["is_raw_unassigned"].mean()),
        "largest_module_fraction": (
            float(proper["assigned_feature_count"].max() / total_assigned)
            if total_assigned
            else None
        ),
        "median_assigned_abs_kME": float(median_kme) if pd.notna(median_kme) else None,
        "collapsed_to_one_proper_module": bool(len(proper) <= 1),
        "single_module_fallback_reason": fallback_reason,
        "dense_adjacency_TOM_persisted": False,
        "model_object_persisted": False,
        "source_matrix": matrix_path.resolve().as_posix(),
        "sample_info": sample_info_path.resolve().as_posix(),
        "runtime_seconds": float(time.perf_counter() - started),
        "findModules_run": True,
    }
    write_json(summary_path, payload)
    output_paths = [
        root / "module_assignments.csv.gz",
        root / "module_summary.csv",
        root / "module_eigengenes.csv",
        root / "dynamic_module_eigengenes.csv",
        root / "module_eigengene_linkage.csv",
        root / "soft_threshold_table.csv",
        root / "gene_tree_linkage.csv.gz",
        sample_info_path,
        summary_path,
        *sorted(figures.glob("*.png")),
    ]
    write_run_manifest(
        manifest_path,
        pipeline_name="transcriptomics_after_outlier_removal_find_modules",
        run_root=run_root,
        parameters={
            "condition": CONDITION,
            "automatic_power": True,
            "deepSplit": 2,
            "minModuleSize": payload["minModuleSize"],
            "MEDissThres": payload["MEDissThres"],
            "pamRespectsDendro": False,
        },
        generation_rules={
            "sample_exclusion": "Ng_MOG1 and Ng_HD3 only",
            "sample_tree_cut": "none; cut=inf preprocess matrix used directly",
            "additional_filtering": "none",
            "expression_transform": "none",
            "power": "PyWGCNA pickSoftThreshold automatic selection per network",
            "storage": "compact tables/linkage/figures; no dense adjacency, TOM, or model pickle",
        },
        input_paths=[matrix_path, preprocess_path, config.clinical_path],
        output_paths=output_paths,
        code_paths=[
            Path(__file__),
            Path(__file__).parent / "run_transcriptomics_variance_find_modules.py",
            Path(__file__).parent / "pipeline_provenance.py",
        ],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )
    result = compact(network, strategy, payload)
    model.adjacency = None
    model.TOM = None
    del model, matrix, metadata, sample_info, assignments, module_summary, raw_dynamic
    plt.close("all")
    gc.collect()
    return result


def collect_summary(config: PipelineConfig) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for network, strategy in ordered_networks(config):
        path = output_root(network) / "analysis_summary.json"
        if path.exists():
            rows.append(compact(network, strategy, json.loads(path.read_text(encoding="utf-8"))))
    return pd.DataFrame(rows)


def write_gallery(config: PipelineConfig, summary: pd.DataFrame) -> Path:
    stage = config.experiment_root / "04_variance_filtered/transcriptomics"
    path = stage / "outlier_removal_find_modules_gallery.html"
    cards: list[str] = []
    for row in summary.itertuples(index=False):
        network = next(
            network
            for network, strategy in ordered_networks(config)
            if network.network_id == row.network_id and strategy == row.strategy
        )
        root = output_root(network)
        image = root / "figures/module_eigengene_dendrogram.png"
        relative = image.relative_to(stage).as_posix() if image.exists() else ""
        cards.append(
            "<article><b>" + html.escape(str(row.cell_type)) + "</b><small>" +
            html.escape(f"{row.specimen} | {row.threshold} | {row.strategy} | {row.mode} | {row.status}") +
            "</small>" +
            (f"<a href='{relative}'><img loading='lazy' src='{relative}'></a>" if relative else "") +
            "</article>"
        )
    path.write_text(
        "<!doctype html><meta charset='utf-8'><title>After-outlier findModules</title>"
        "<style>body{font:14px Arial;margin:18px;color:#17233a}main{display:grid;grid-template-columns:repeat(auto-fill,minmax(420px,1fr));gap:12px}article{border:1px solid #d7dde7;padding:10px}small{display:block;color:#64748b;margin:5px 0}img{width:100%;height:330px;object-fit:contain}</style>"
        "<h1>findModules after Ng_MOG1 and Ng_HD3 removal</h1><p>No finite sample-tree cut was applied.</p><main>" + "".join(cards) + "</main>",
        encoding="utf-8",
    )
    return path


def run(args: argparse.Namespace) -> None:
    run_root = args.run_root.resolve()
    config = PipelineConfig(run_root=run_root)
    entries = ordered_networks(config)
    if args.only_network:
        entries = [item for item in entries if item[0].network_id == args.only_network]
    if args.limit is not None:
        entries = entries[: args.limit]
    if args.gallery_only:
        summary = collect_summary(config)
        summary.to_csv(
            config.experiment_root / "04_variance_filtered/transcriptomics/outlier_removal_find_modules_summary.csv",
            index=False,
        )
        print(f"Updated {write_gallery(config, summary)}", flush=True)
        return
    if args.worker_count < 1 or not 0 <= args.worker_index < args.worker_count:
        raise ValueError("worker index/count are invalid")
    entries = entries[args.worker_index :: args.worker_count]
    clinical = pd.read_csv(config.clinical_path)
    cohort_map = build_source_cohort_map(config.run_root)
    for position, (network, strategy) in enumerate(entries, start=1):
        print(
            f"[after-outlier findModules {position}/{len(entries)}] {strategy} / {network.network_id}",
            flush=True,
        )
        try:
            result = process_network(run_root, config, network, strategy, clinical, cohort_map)
            print(
                f"  {result['status']} | modules={result.get('final_module_count')} | "
                f"power={result.get('selected_power')}",
                flush=True,
            )
        except BaseException as error:
            root = output_root(network)
            payload = {
                "status": "failed_findModules",
                "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                "error_type": type(error).__name__,
                "error": str(error),
                "findModules_run": True,
                "finite_sample_cut_applied": False,
            }
            write_json(root / "analysis_summary.json", payload)
            print(f"  FAILED: {type(error).__name__}: {error}", flush=True)
        gc.collect()
    print(
        f"Completed worker {args.worker_index + 1}/{args.worker_count}; run --gallery-only after all workers.",
        flush=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    parser.add_argument("--only-network")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--gallery-only", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
