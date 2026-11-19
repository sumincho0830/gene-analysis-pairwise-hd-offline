"""Recalculate selected transcriptomics ME merges at MEDissThres=0.10.

Only the module-eigengene merge parameter changes relative to the saved
baseline.  Therefore this runner reuses the saved dynamic tree-cut modules and
recalculates the official PyWGCNA ``mergeCloseModules`` stage, rather than
recomputing the unchanged adjacency, TOM, gene tree, and dynamic tree cut.

For every selected network, the same merge-stage reconstruction is first run
at MEDissThres=0.20.  The recalculation proceeds only if that partition exactly
matches the saved baseline partition (module color names may be relabelled, but
the dynamic-module grouping must be identical).
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
import time
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_provenance import write_run_manifest
from run_transcriptomics_variance_find_modules import (
    compact_summary,
    export_tables,
    finish_module_merge_after_pywgcna_exit,
    generate_dynamic_me_dendrogram,
    load_model,
    network_root,
    plot_gene_dendrogram_modules,
    plot_module_sizes,
    plot_soft_threshold,
    validate_input_model,
    write_json,
)


BASELINE_CONDITION = "01_default_ds2_mms50_med0p20"
OUTPUT_CONDITION = "02_predicted_ds2_mms50_med0p10"
BASELINE_THRESHOLD = 0.20
CANDIDATE_THRESHOLD = 0.10


def normalized_text(value: Any) -> str:
    if pd.isna(value) or str(value).strip() == "":
        return "cell_integrated"
    return str(value)


def network_key(row: pd.Series) -> str:
    return "|".join(
        normalized_text(row[column])
        for column in ("threshold", "strategy", "mode", "specimen", "cell_type")
    )


def partition_signature(
    dynamic_colors: np.ndarray, final_colors: np.ndarray
) -> tuple[tuple[str, ...], ...]:
    groups: dict[str, set[str]] = {}
    for dynamic, final in zip(dynamic_colors, final_colors, strict=True):
        groups.setdefault(str(final), set()).add(str(dynamic))
    return tuple(sorted(tuple(sorted(members)) for members in groups.values()))


def set_model_merge_result(
    model: Any,
    dynamic_colors: np.ndarray,
    threshold: float,
) -> str | None:
    """Run official PyWGCNA ME merging and attach the result to ``model``."""

    from PyWGCNA.wgcna import WGCNA

    expression = model.datExpr.to_df()
    dynamic = np.asarray(dynamic_colors, dtype=object).copy()
    model.datExpr.var["dynamicColors"] = dynamic.astype(str)
    model.MEDissThres = float(threshold)
    fallback_reason: str | None = None
    try:
        merge = WGCNA.mergeCloseModules(
            expression,
            dynamic,
            cutHeight=float(threshold),
        )
        if not bool(merge.get("allOK", False)):
            raise RuntimeError("PyWGCNA mergeCloseModules returned allOK=False")
        merged = np.asarray(merge["colors"], dtype=str)
        model.datExpr.var["moduleColors"] = merged
        color_order = np.unique(merged).tolist()
        model.datExpr.var["moduleLabels"] = [
            color_order.index(value) for value in merged
        ]
        model.datME = WGCNA.moduleEigengenes(
            expression,
            model.datExpr.var["moduleColors"],
        )["eigengenes"]
        if "MEgrey" in model.datME.columns:
            model.datME = model.datME.drop(columns=["MEgrey"])
        model.MEs = (
            WGCNA.orderMEs(model.datME)
            if model.datME.shape[1] > 1
            else model.datME.copy()
        )
    except (Exception, SystemExit) as error:
        fallback_reason = f"{type(error).__name__}: {error}"
        model.datExpr.var["dynamicColors"] = np.asarray(
            dynamic_colors, dtype=str
        )
        model.MEDissThres = float(threshold)
        finish_module_merge_after_pywgcna_exit(model)
    return fallback_reason


def prepare_model_from_baseline(
    model: Any,
    baseline_root: Path,
    baseline_assignments: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame, dict[str, Any]]:
    baseline_summary = json.loads(
        (baseline_root / "analysis_summary.json").read_text(encoding="utf-8")
    )
    expression = model.datExpr.to_df()
    expression.columns = expression.columns.astype(str)
    assignments = baseline_assignments.copy()
    assignments["gene_id"] = assignments["gene_id"].astype(str)
    if assignments["gene_id"].duplicated().any():
        raise ValueError("Duplicate gene_id in baseline module assignments")
    assignments = assignments.set_index("gene_id").reindex(expression.columns)
    if assignments.isna().all(axis=1).any():
        missing = assignments.index[assignments.isna().all(axis=1)].tolist()
        raise ValueError(f"Missing {len(missing)} expression genes in assignments")
    assignments.index.name = "gene_id"
    assignments = assignments.reset_index()
    if assignments["gene_id"].tolist() != expression.columns.tolist():
        raise ValueError("Baseline assignments could not be aligned to expression genes")

    dynamic_colors = assignments["dynamic_color"].astype(str).to_numpy()
    raw_dynamic = assignments[
        ["raw_dynamic_label", "ordered_dynamic_label"]
    ].rename(
        columns={"raw_dynamic_label": "Name", "ordered_dynamic_label": "Value"}
    )
    raw_dynamic = raw_dynamic.astype(int)

    model.power = int(baseline_summary["selected_power"])
    model.sft = pd.read_csv(baseline_root / "soft_threshold_table.csv")
    model.geneTree = pd.read_csv(
        baseline_root / "gene_tree_linkage.csv.gz"
    ).to_numpy(dtype=float)
    model.deepSplit = int(baseline_summary["deepSplit"])
    model.minModuleSize = int(baseline_summary["minModuleSize"])
    model.save = False
    return dynamic_colors, raw_dynamic, baseline_summary


def process_network(
    run_root: Path,
    stage_root: Path,
    row: pd.Series,
    predicted_row: pd.Series,
    force: bool,
) -> dict[str, Any]:
    root = network_root(stage_root, row)
    baseline_root = root / "05_find_modules" / BASELINE_CONDITION
    output_root = root / "05_find_modules" / OUTPUT_CONDITION
    figures_root = output_root / "figures"
    summary_path = output_root / "analysis_summary.json"
    manifest_path = output_root / "run_manifest.json"
    if not force and summary_path.exists() and manifest_path.exists():
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        if payload.get("status") == "complete_findModules_merge_recalculation":
            return compact_summary(row, payload)

    preprocess_summary_path = root / "03_preprocess/02_selected_cut/preprocess_summary.json"
    preprocess = json.loads(preprocess_summary_path.read_text(encoding="utf-8"))
    if preprocess.get("status") != "complete_selected_cut":
        raise ValueError(f"Selected-cut preprocess is not complete: {preprocess.get('status')}")
    object_path = Path(str(preprocess["model_object"]))
    baseline_assignments_path = baseline_root / "module_assignments.csv.gz"
    baseline_assignments = pd.read_csv(baseline_assignments_path)

    output_root.mkdir(parents=True, exist_ok=True)
    figures_root.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    model = load_model(object_path)
    validate_input_model(
        model,
        expected_samples=int(preprocess["samples_after_preprocess"]),
        expected_genes=int(preprocess["genes_after_preprocess"]),
    )
    dynamic_colors, raw_dynamic, baseline_summary = prepare_model_from_baseline(
        model,
        baseline_root,
        baseline_assignments,
    )

    # Gate: reconstruct the saved baseline merge before changing its threshold.
    baseline_fallback = set_model_merge_result(
        model,
        dynamic_colors,
        BASELINE_THRESHOLD,
    )
    reconstructed_baseline = model.datExpr.var["moduleColors"].astype(str).to_numpy()
    saved_baseline = baseline_assignments.set_index("gene_id").loc[
        model.datExpr.to_df().columns.astype(str), "module_color"
    ].astype(str).to_numpy()
    baseline_partition_match = partition_signature(
        dynamic_colors, reconstructed_baseline
    ) == partition_signature(dynamic_colors, saved_baseline)
    baseline_color_match = bool(np.array_equal(reconstructed_baseline, saved_baseline))
    if not baseline_partition_match:
        raise RuntimeError(
            "Official merge-stage reconstruction at MEDissThres=0.20 did not "
            "reproduce the saved baseline partition"
        )

    candidate_fallback = set_model_merge_result(
        model,
        dynamic_colors,
        CANDIDATE_THRESHOLD,
    )
    assignments, module_summary = export_tables(model, raw_dynamic, output_root)
    generate_dynamic_me_dendrogram(
        model.datExpr.to_df(),
        assignments,
        output_root,
        me_diss_threshold=CANDIDATE_THRESHOLD,
    )

    label = (
        f"{row['threshold']} / {row['strategy']} / {row['mode']} / "
        f"{row['specimen']} / {row['cell_type']}"
    )
    plot_gene_dendrogram_modules(
        model,
        figures_root / "gene_dendrogram_modules.png",
        f"Transcriptomics modules at MEDissThres=0.10: {label}",
    )
    plot_soft_threshold(
        model,
        figures_root / "soft_threshold.png",
        f"Transcriptomics soft-threshold (reused baseline): {label}",
    )
    plot_module_sizes(
        module_summary,
        figures_root / "module_sizes.png",
        f"Transcriptomics module sizes at MEDissThres=0.10: {label}",
    )

    proper = module_summary.loc[~module_summary["is_fully_raw_unassigned"]]
    total_assigned = int(proper["assigned_feature_count"].sum())
    largest_fraction = (
        float(proper["assigned_feature_count"].max() / total_assigned)
        if total_assigned
        else None
    )
    median_kme = assignments.loc[
        ~assignments["is_raw_unassigned"].astype(bool), "assigned_abs_kME"
    ].median()
    actual_final_count = int(len(proper))
    current_final_count = int(baseline_summary["final_module_count"])
    payload = {
        "status": "complete_findModules_merge_recalculation",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "condition": OUTPUT_CONDITION,
        "calculation_scope": "official_PyWGCNA_mergeCloseModules_recalculation",
        "samples": int(model.datExpr.shape[0]),
        "genes": int(model.datExpr.shape[1]),
        "feature_filtering_before_findModules": "none",
        "expression_transform_before_findModules": "none",
        "selected_power": int(model.power),
        "selected_power_source": "reused_identical_baseline_network",
        "networkType": str(model.networkType),
        "TOMType": str(model.TOMType),
        "deepSplit": 2,
        "minModuleSize": 50,
        "MEDissThres": CANDIDATE_THRESHOLD,
        "pamRespectsDendro": False,
        "dynamic_module_count": int(
            assignments.loc[
                ~assignments["is_raw_unassigned"].astype(bool), "raw_dynamic_label"
            ].nunique()
        ),
        "baseline_final_module_count": current_final_count,
        "predicted_diagnostic_module_count": int(
            predicted_row["diagnostic_candidate_module_count"]
        ),
        "actual_final_module_count": actual_final_count,
        "actual_additional_modules": actual_final_count - current_final_count,
        "prediction_minus_actual": int(
            predicted_row["diagnostic_candidate_module_count"]
        ) - actual_final_count,
        "final_module_count": actual_final_count,
        "raw_unassigned_gene_count": int(assignments["is_raw_unassigned"].sum()),
        "raw_unassigned_fraction": float(assignments["is_raw_unassigned"].mean()),
        "largest_module_fraction": largest_fraction,
        "median_assigned_abs_kME": (
            float(median_kme) if pd.notna(median_kme) else None
        ),
        "collapsed_to_one_proper_module": bool(actual_final_count <= 1),
        "baseline_partition_reproduced": bool(baseline_partition_match),
        "baseline_color_labels_reproduced": bool(baseline_color_match),
        "baseline_merge_fallback_reason": baseline_fallback,
        "candidate_merge_fallback_reason": candidate_fallback,
        "upstream_network_recomputed": False,
        "upstream_network_reuse_reason": (
            "deepSplit, minModuleSize, power, adjacency, TOM, gene tree, and "
            "dynamic modules are unchanged; only MEDissThres is recalculated"
        ),
        "dense_adjacency_TOM_persisted": False,
        "model_object_persisted": False,
        "runtime_seconds": float(time.perf_counter() - started),
        "findModules_run": True,
    }
    write_json(summary_path, payload)
    output_paths = [
        output_root / "module_assignments.csv.gz",
        output_root / "module_summary.csv",
        output_root / "module_eigengenes.csv",
        output_root / "dynamic_module_eigengenes.csv",
        output_root / "module_eigengene_linkage.csv",
        output_root / "soft_threshold_table.csv",
        output_root / "gene_tree_linkage.csv.gz",
        summary_path,
        *sorted(figures_root.glob("*.png")),
    ]
    write_run_manifest(
        manifest_path,
        pipeline_name="transcriptomics_predicted_MEDissThres_merge_recalculation",
        run_root=run_root,
        parameters={
            "condition": OUTPUT_CONDITION,
            "deepSplit": 2,
            "minModuleSize": 50,
            "MEDissThres": CANDIDATE_THRESHOLD,
            "pamRespectsDendro": False,
            "selected_power": int(model.power),
            "power_source": "saved baseline",
        },
        generation_rules={
            "selection": "networks with a saved 0.10 ME merge candidate",
            "input": "selected-cut expression plus saved baseline dynamic modules",
            "calculation": "official PyWGCNA mergeCloseModules at cutHeight=0.10",
            "baseline_gate": (
                "same merge-stage reconstruction at 0.20 must reproduce the "
                "saved baseline dynamic-module partition"
            ),
            "unchanged_upstream": (
                "adjacency, TOM, gene tree, and dynamic tree cut are reused "
                "because their parameters are unchanged"
            ),
            "additional_filtering_or_transform": "none",
            "storage": "same compact tables and figures as baseline findModules",
        },
        input_paths=[
            object_path,
            preprocess_summary_path,
            baseline_assignments_path,
            baseline_root / "analysis_summary.json",
            baseline_root / "soft_threshold_table.csv",
            baseline_root / "gene_tree_linkage.csv.gz",
        ],
        output_paths=output_paths,
        code_paths=[
            Path(__file__),
            Path(__file__).parent / "run_transcriptomics_variance_find_modules.py",
            Path(__file__).parent / "pipeline_provenance.py",
        ],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )
    result = compact_summary(row, payload)
    model.adjacency = None
    model.TOM = None
    del model, assignments, module_summary, raw_dynamic, baseline_assignments
    plt.close("all")
    gc.collect()
    return result


def run(
    run_root: Path,
    force: bool,
    limit: int | None,
    worker_index: int,
    worker_count: int,
) -> None:
    stage_root = run_root / "experiments/04_variance_filtered/transcriptomics"
    diagnostic_root = stage_root / "07_module_merge_diagnostics"
    candidates_path = diagnostic_root / "network_candidates.csv"
    candidates = pd.read_csv(candidates_path)
    candidates = candidates.loc[candidates["potential_additional_modules"] > 0].copy()
    if limit is not None:
        candidates = candidates.head(limit).copy()
    if worker_count < 1:
        raise ValueError("worker_count must be at least 1")
    if not 0 <= worker_index < worker_count:
        raise ValueError("worker_index must satisfy 0 <= worker_index < worker_count")
    candidates = candidates.iloc[worker_index::worker_count].copy()

    registry = pd.read_csv(stage_root / "find_modules_summary.csv")
    registry["_key"] = registry.apply(network_key, axis=1)
    registry = registry.set_index("_key", drop=False)
    results: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for position, (_, predicted) in enumerate(candidates.iterrows(), start=1):
        key = str(predicted["network_key"])
        print(f"[{position}/{len(candidates)}] {key}", flush=True)
        if key not in registry.index:
            errors.append({"network_key": key, "error": "missing_from_registry"})
            continue
        row = registry.loc[key]
        if isinstance(row, pd.DataFrame):
            raise ValueError(f"Duplicate registry key: {key}")
        try:
            result = process_network(
                run_root,
                stage_root,
                row,
                predicted,
                force=force,
            )
            results.append(result)
            print(
                "  actual final modules: "
                f"{result['actual_final_module_count']} "
                f"(baseline {result['baseline_final_module_count']}, "
                f"predicted {result['predicted_diagnostic_module_count']})",
                flush=True,
            )
        except Exception as error:
            errors.append(
                {
                    "network_key": key,
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            root = network_root(stage_root, row) / "05_find_modules" / OUTPUT_CONDITION
            write_json(
                root / "analysis_summary.json",
                {
                    "status": "failed_merge_recalculation",
                    "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                    "error_type": type(error).__name__,
                    "error": str(error),
                },
            )
            print(f"  FAILED: {type(error).__name__}: {error}", flush=True)
        gc.collect()

    result_table = pd.DataFrame(results)
    summary_csv = diagnostic_root / "actual_rerun_summary.csv"
    errors_path = diagnostic_root / "actual_rerun_errors.json"
    analysis_path = diagnostic_root / "actual_rerun_analysis_summary.json"
    manifest_path = diagnostic_root / "actual_rerun_manifest.json"
    result_table.to_csv(summary_csv, index=False)
    write_json(errors_path, {"errors": errors})
    complete = result_table.loc[
        result_table.get("status", pd.Series(dtype=str))
        == "complete_findModules_merge_recalculation"
    ]
    summary = {
        "status": "complete" if not errors else "complete_with_errors",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selected_networks": int(len(candidates)),
        "completed_networks": int(len(complete)),
        "failed_networks": int(len(errors)),
        "baseline_partition_validation_passed": int(
            complete.get("baseline_partition_reproduced", pd.Series(dtype=bool))
            .astype(bool)
            .sum()
        ),
        "baseline_color_label_matches": int(
            complete.get("baseline_color_labels_reproduced", pd.Series(dtype=bool))
            .astype(bool)
            .sum()
        ),
        "collapsed_before": int(
            (complete.get("baseline_final_module_count", pd.Series(dtype=float)) <= 1).sum()
        ),
        "collapsed_after": int(
            (complete.get("actual_final_module_count", pd.Series(dtype=float)) <= 1).sum()
        ),
        "total_baseline_modules": int(
            complete.get("baseline_final_module_count", pd.Series(dtype=float)).sum()
        ),
        "total_predicted_diagnostic_modules": int(
            complete.get(
                "predicted_diagnostic_module_count", pd.Series(dtype=float)
            ).sum()
        ),
        "total_actual_modules": int(
            complete.get("actual_final_module_count", pd.Series(dtype=float)).sum()
        ),
        "median_actual_additional_modules": (
            float(complete["actual_additional_modules"].median())
            if not complete.empty
            else None
        ),
        "calculation_scope": (
            "official PyWGCNA mergeCloseModules recalculation from identical "
            "saved dynamic modules; no upstream recomputation"
        ),
    }
    write_json(analysis_path, summary)
    write_run_manifest(
        manifest_path,
        pipeline_name="all_transcriptomics_predicted_ME_merge_recalculations",
        run_root=run_root,
        parameters={
            "condition": OUTPUT_CONDITION,
            "baseline_MEDissThres": BASELINE_THRESHOLD,
            "candidate_MEDissThres": CANDIDATE_THRESHOLD,
            "force": force,
            "limit": limit,
        },
        generation_rules={
            "selection": "network_candidates.csv potential_additional_modules > 0",
            "baseline_validation": "required before every candidate merge",
            "calculation": "official mergeCloseModules only; unchanged upstream reused",
        },
        input_paths=[
            candidates_path,
            stage_root / "find_modules_summary.csv",
        ],
        output_paths=[summary_csv, errors_path, analysis_path],
        code_paths=[
            Path(__file__),
            Path(__file__).parent / "run_transcriptomics_variance_find_modules.py",
            Path(__file__).parent / "pipeline_provenance.py",
        ],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(
        arguments.run_root.resolve(),
        force=arguments.force,
        limit=arguments.limit,
        worker_index=arguments.worker_index,
        worker_count=arguments.worker_count,
    )
