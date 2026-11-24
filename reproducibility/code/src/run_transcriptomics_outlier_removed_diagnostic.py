"""Re-run transcriptomics QC and cut=inf preprocess after two exclusions.

The existing variance-filtered results remain the immutable
``Before Outlier Removal`` branch.  This script removes the two confirmed
biological subjects, regenerates the ten descriptive/distribution outputs,
and runs PyWGCNA ``preprocess(cut=inf)`` only.  A sample dendrogram is retained
as diagnostic evidence, but no selected-cut stage is created.  The cut=inf
matrix proceeds directly to findModules.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform

from pipeline_provenance import write_run_manifest
from run_transcriptomics_stage1 import (
    PipelineConfig,
    grouped_distance_robust_z,
    pywgcna_pandas_compatibility,
)
from run_transcriptomics_variance_filtered import (
    DIAGNOSTIC_ROBUST_Z,
    STRATEGIES,
    Network,
    assert_matrix_equal,
    audit_pre_cut,
    link_or_write_matrix_gz,
    load_networks,
    network_metadata,
    safe_write_json,
    write_matrix_gz,
)
from transcriptomics_qc import (
    build_source_cohort_map,
    read_expression_matrix,
    write_qc_outputs,
)


EXCLUDED_SUBJECTS = ("Ng_MOG1", "Ng_HD3")
CAUTION_REMOVAL_FRACTION = 0.05
SUBJECT_PATTERNS = {
    subject: re.compile(rf"(?:^|_){re.escape(subject)}(?:_|$)")
    for subject in EXCLUDED_SUBJECTS
}


def source_sample_id(row_id: str) -> str:
    """Return the biological sample ID from separate or integrated row IDs."""
    return str(row_id).split("__CELL__", 1)[0]


def excluded_subject(row_id: str) -> str | None:
    sample_id = source_sample_id(row_id)
    for subject, pattern in SUBJECT_PATTERNS.items():
        if pattern.search(sample_id):
            return subject
    return None


def exclusion_table(matrix: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row_id in matrix.index.astype(str):
        subject = excluded_subject(row_id)
        if subject is None:
            continue
        meta = metadata.loc[row_id]
        rows.append(
            {
                "matrix_row_id": row_id,
                "source_sample_id": source_sample_id(row_id),
                "excluded_subject": subject,
                "specimen": meta.get("specimen"),
                "cell_type": meta.get("cell_type"),
                "cohort": meta.get("cohort"),
                "diagnosis": meta.get("diagnosis"),
                "reason": "confirmed severe PCA outlier; manual decision",
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "matrix_row_id",
            "source_sample_id",
            "excluded_subject",
            "specimen",
            "cell_type",
            "cohort",
            "diagnosis",
            "reason",
        ],
    )


def linkage_gap_table(
    linkage: np.ndarray,
    sample_ids: list[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Record every merge-height gap and suggest the global largest midpoint.

    The suggestion follows the deliberately simple visual rule chosen for the
    manual review: take the largest difference between consecutive linkage
    heights and place the line at its midpoint.  A warning is attached when
    retaining the unique largest cluster would remove more than five percent
    of rows, or when the largest cluster is tied.  The suggestion is never
    applied automatically.
    """
    heights = linkage[:, 2]
    gaps = np.diff(heights)
    rows: list[dict[str, Any]] = []
    max_removed = max(1, int(len(sample_ids) * CAUTION_REMOVAL_FRACTION))
    for index, gap in enumerate(gaps):
        lower = float(heights[index])
        upper = float(heights[index + 1])
        cut = (lower + upper) / 2.0
        labels = fcluster(linkage, t=cut, criterion="distance")
        unique, counts = np.unique(labels, return_counts=True)
        largest = int(counts.max())
        unique_largest = int(np.sum(counts == largest)) == 1
        main_label = int(unique[int(np.argmax(counts))]) if unique_largest else None
        removed = (
            sorted(
                sample_ids[position]
                for position, label in enumerate(labels)
                if int(label) != main_label
            )
            if main_label is not None
            else list(sample_ids)
        )
        relative_gap = float(gap / upper) if upper > 0 else 0.0
        caution = (not unique_largest) or len(removed) > max_removed
        rows.append(
            {
                "merge_index_lower": int(index),
                "lower_height": lower,
                "upper_height": upper,
                "cut_midpoint": cut,
                "absolute_gap": float(gap),
                "relative_gap": relative_gap,
                "unique_largest_cluster": bool(unique_largest),
                "removed_count_if_largest_cluster_retained": len(removed),
                "max_removed_for_review": max_removed,
                "requires_caution": caution,
                "removed_sample_ids_if_largest_cluster_retained": "; ".join(removed),
            }
        )
    table = pd.DataFrame(rows)
    if table.empty or not np.isfinite(table["absolute_gap"]).any():
        return table, {
            "status": "no_finite_linkage_gap",
            "suggested_cut": None,
            "applied": False,
        }
    best = table.sort_values(
        ["absolute_gap", "relative_gap"],
        ascending=[False, False],
        kind="mergesort",
    ).iloc[0]
    table["selected_review_suggestion"] = False
    table.loc[best.name, "selected_review_suggestion"] = True
    removed = [
        value.strip()
        for value in str(
            best["removed_sample_ids_if_largest_cluster_retained"]
        ).split(";")
        if value.strip()
    ]
    return table, {
        "status": "review_required",
        "rule": "midpoint of the global largest consecutive linkage-height gap",
        "caution_removal_fraction": CAUTION_REMOVAL_FRACTION,
        "suggested_cut": float(best["cut_midpoint"]),
        "lower_height": float(best["lower_height"]),
        "upper_height": float(best["upper_height"]),
        "absolute_gap": float(best["absolute_gap"]),
        "relative_gap": float(best["relative_gap"]),
        "removed_count_if_largest_cluster_retained": int(
            best["removed_count_if_largest_cluster_retained"]
        ),
        "removed_sample_ids_if_largest_cluster_retained": removed,
        "unique_largest_cluster": bool(best["unique_largest_cluster"]),
        "requires_caution": bool(best["requires_caution"]),
        "applied": False,
        "decision": "pending_manual_review",
    }


def plot_diagnostic(
    linkage: np.ndarray,
    labels: list[str],
    recommendation: dict[str, Any],
    output_path: Path,
    title: str,
) -> list[int]:
    width = max(18.0, min(50.0, len(labels) / 8.0))
    figure, axis = plt.subplots(figsize=(width, 9), facecolor="white")
    suggested = recommendation.get("suggested_cut")
    result = dendrogram(
        linkage,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=5,
        color_threshold=0,
        above_threshold_color="#4C78A8",
        ax=axis,
    )
    if suggested is not None and np.isfinite(float(suggested)):
        axis.axhline(
            float(suggested),
            color="#D62728",
            linewidth=1.4,
            linestyle="--",
            label=f"largest-gap suggestion = {float(suggested):.6g}",
        )
        axis.legend(loc="upper right", frameon=True)
        suffix = f"largest-gap suggestion={float(suggested):.6g}; not applied"
    else:
        suffix = "no finite suggestion; no cut applied"
    axis.set_title(f"{title}\n{suffix}")
    axis.set_xlabel("Samples")
    axis.set_ylabel("Euclidean distance; average linkage")
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=210, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [int(value) for value in result["leaves"]]


def write_pending_selected_cut_page(
    output_path: Path,
    network: Network,
    strategy: str,
    recommendation: dict[str, Any],
) -> None:
    suggested = recommendation.get("suggested_cut")
    suggested_text = "None" if suggested is None else f"{float(suggested):.6g}"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "<!doctype html><meta charset='utf-8'><title>Selected cut pending</title>"
        "<style>body{font:16px/1.55 Arial,sans-serif;max-width:820px;margin:48px auto;"
        "padding:0 24px;color:#17233a}code{background:#f2f4f7;padding:2px 5px}</style>"
        f"<h1>Selected cut pending manual review</h1>"
        f"<p><b>Network:</b> {network.network_id}</p>"
        f"<p><b>Strategy:</b> {strategy}</p>"
        f"<p><b>Largest-gap suggestion:</b> <code>{suggested_text}</code></p>"
        "<p>The suggestion is diagnostic evidence only. No finite cut has been "
        "applied and no selected-cut matrix has been generated.</p>",
        encoding="utf-8",
    )


def run_diagnostic(
    network: Network,
    strategy: str,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    selected_matrix_path: Path,
) -> tuple[dict[str, Any], list[Path]]:
    from PyWGCNA import WGCNA

    branch_root = network.output_root / "outlier_removal"
    diagnostic_root = branch_root / "03_preprocess/01_cut_diagnostic"
    diagnostic_root.mkdir(parents=True, exist_ok=True)

    if matrix.shape[0] < 4:
        summary = {
            "status": "not_run_insufficient_samples",
            "input_samples": int(matrix.shape[0]),
            "input_genes": int(matrix.shape[1]),
            "cut_applied": False,
        }
        safe_write_json(diagnostic_root / "preprocess_summary.json", summary)
        return summary, [diagnostic_root / "preprocess_summary.json"]

    sample_info = metadata.loc[
        matrix.index,
        [
            "project_id",
            "patient_id",
            "cohort",
            "specimen",
            "cell_type",
            "diagnosis",
            "sex",
            "age",
        ],
    ].copy()
    model = WGCNA(
        name=f"TRANSCRIPTOMICS_{network.network_id}_{strategy}_after_outlier_removal_cut_inf",
        geneExp=matrix,
        sampleInfo=sample_info,
        species="Human",
        save=False,
        outputPath=diagnostic_root.resolve().as_posix() + "/",
        figureType="png",
        cut=float("inf"),
    )
    expected, tpm_keep, good_genes, good_samples = audit_pre_cut(matrix, model)
    with pywgcna_pandas_compatibility():
        model.preprocess(show=False)
    before_cut = model.datExpr.to_df().copy()
    assert_matrix_equal(
        before_cut,
        expected,
        f"{network.network_id}/{strategy}/after_outlier_removal/cut_inf",
    )
    before_path = diagnostic_root / "matrix_before_cut.csv.gz"
    link_or_write_matrix_gz(
        before_cut,
        before_path,
        link_source=selected_matrix_path,
        link_source_matrix=matrix,
    )

    removed_features = [
        {
            "gene_symbol": gene,
            "removal_stage": "TPMcutoff",
            "reason": f"no sample value greater than TPMcutoff={model.TPMcutoff}",
        }
        for gene in matrix.columns[~tpm_keep]
    ]
    after_tpm = matrix.loc[:, tpm_keep]
    removed_features.extend(
        {
            "gene_symbol": gene,
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        }
        for gene in after_tpm.columns[~good_genes]
    )
    removed_features_path = diagnostic_root / "removed_features.csv"
    pd.DataFrame(
        removed_features,
        columns=["gene_symbol", "removal_stage", "reason"],
    ).to_csv(removed_features_path, index=False)
    removed_samples_path = diagnostic_root / "removed_samples.csv"
    pd.DataFrame(
        {
            "sample_id": after_tpm.index[~good_samples].astype(str),
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        }
    ).to_csv(removed_samples_path, index=False)

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
        metrics["diagnostic_outlier_candidate"], "review_only", "none"
    )
    metrics_path = diagnostic_root / "sample_distance_metrics.csv"
    metrics.reset_index(drop=True).to_csv(metrics_path, index=False)
    flagged_path = diagnostic_root / "outlier_candidates_before_cut.csv"
    metrics.loc[metrics["diagnostic_outlier_candidate"]].reset_index(
        drop=True
    ).to_csv(flagged_path, index=False)

    candidate_table, recommendation = linkage_gap_table(
        linkage, before_cut.index.astype(str).tolist()
    )
    candidates_path = diagnostic_root / "cut_candidates.csv"
    candidate_table.to_csv(candidates_path, index=False)
    recommendation_path = diagnostic_root / "cut_recommendation.json"
    safe_write_json(recommendation_path, recommendation)
    dendrogram_path = diagnostic_root / "dendrogram_cut_diagnostic.png"
    leaf_order = plot_diagnostic(
        linkage,
        before_cut.index.astype(str).tolist(),
        recommendation,
        dendrogram_path,
        f"After outlier removal: {network.network_id} / {strategy}",
    )
    leaf_path = diagnostic_root / "dendrogram_leaf_order.csv"
    pd.DataFrame(
        {
            "leaf_position_1based": np.arange(1, len(leaf_order) + 1),
            "sample_id": before_cut.index[leaf_order].astype(str),
        }
    ).to_csv(leaf_path, index=False)

    summary = {
        "status": "complete_preprocess_cut_inf",
        "network_id": network.network_id,
        "strategy": strategy,
        "input_samples_after_outlier_removal": int(matrix.shape[0]),
        "input_genes": int(matrix.shape[1]),
        "TPMcutoff": float(model.TPMcutoff),
        "genes_removed_by_TPMcutoff": int((~tpm_keep).sum()),
        "genes_removed_by_goodSamplesGenes": int((~good_genes).sum()),
        "samples_removed_by_goodSamplesGenes": int((~good_samples).sum()),
        "matrix_before_cut_samples": int(before_cut.shape[0]),
        "matrix_before_cut_genes": int(before_cut.shape[1]),
        "diagnostic_outlier_candidates": int(
            metrics["diagnostic_outlier_candidate"].sum()
        ),
        "largest_gap_suggested_cut": recommendation.get("suggested_cut"),
        "finite_cut_applied": False,
        "downstream_input_status": "approved_without_finite_sample_cut",
        "diagnostic_model_object_saved": False,
    }
    summary_path = diagnostic_root / "preprocess_summary.json"
    safe_write_json(summary_path, summary)
    del model
    gc.collect()
    return summary, [
        before_path,
        removed_features_path,
        removed_samples_path,
        linkage_path,
        metrics_path,
        flagged_path,
        candidates_path,
        recommendation_path,
        dendrogram_path,
        leaf_path,
        summary_path,
    ]


def run(
    config: PipelineConfig,
    *,
    strategies: tuple[str, ...],
    thresholds: tuple[str, ...],
    only_network: str | None,
    resume: bool,
) -> None:
    networks = [
        network
        for network in load_networks(config, strategies)
        if network.threshold in thresholds
        and (only_network is None or network.network_id == only_network)
    ]
    if not networks:
        raise ValueError("No networks matched the requested filters")

    clinical = pd.read_csv(config.clinical_path)
    cohort_map = build_source_cohort_map(config.run_root)
    summary_rows: list[dict[str, Any]] = []
    for position, network in enumerate(networks, start=1):
        strategy_matches = [value for value in STRATEGIES if value in network.output_root.parts]
        if len(strategy_matches) != 1:
            raise ValueError(f"Cannot resolve strategy from {network.output_root}")
        strategy = strategy_matches[0]
        branch_root = network.output_root / "outlier_removal"
        completion_path = branch_root / "03_preprocess/01_cut_diagnostic/preprocess_summary.json"
        if resume and completion_path.exists():
            summary = json.loads(completion_path.read_text(encoding="utf-8"))
            selection_summary_path = branch_root / "01_sample_selection/selection_summary.json"
            selection = (
                json.loads(selection_summary_path.read_text(encoding="utf-8"))
                if selection_summary_path.exists()
                else {}
            )
            summary_rows.append(
                {
                    "threshold": network.threshold,
                    "strategy": strategy,
                    "mode": network.mode,
                    "specimen": network.specimen,
                    "cell_type": network.cell_type,
                    "network_id": network.network_id,
                    "status": summary.get("status"),
                    "input_rows_before_removal": selection.get("input_rows"),
                    "rows_removed": selection.get("rows_removed"),
                    "rows_after_removal": selection.get("rows_after_removal"),
                    "largest_gap_suggested_cut": summary.get("largest_gap_suggested_cut"),
                    "downstream_input_status": summary.get("downstream_input_status"),
                }
            )
            continue

        print(
            f"[outlier-removal {position}/{len(networks)}] "
            f"{strategy} / {network.network_id}",
            flush=True,
        )
        selected_matrix_path = (
            network.output_root
            / "01_feature_selection/selected_expression_matrix.csv.gz"
        )
        if not selected_matrix_path.exists():
            raise FileNotFoundError(selected_matrix_path)
        matrix = read_expression_matrix(selected_matrix_path)
        metadata = network_metadata(network, matrix, clinical, cohort_map)
        excluded = exclusion_table(matrix, metadata)
        remove_mask = matrix.index.astype(str).map(excluded_subject).notna()
        retained = matrix.loc[~remove_mask].copy()
        retained_metadata = metadata.loc[retained.index].copy()
        if any(excluded_subject(value) for value in retained.index.astype(str)):
            raise AssertionError("Confirmed outlier rows remained after filtering")

        selection_root = branch_root / "01_sample_selection"
        selection_root.mkdir(parents=True, exist_ok=True)
        excluded_path = selection_root / "excluded_samples.csv"
        excluded.to_csv(excluded_path, index=False)
        retained_path = selection_root / "retained_samples.csv"
        retained_metadata.reset_index(drop=True).to_csv(retained_path, index=False)
        retained_matrix_path = selection_root / "matrix_after_outlier_removal.csv.gz"
        write_matrix_gz(retained, retained_matrix_path)
        selection_summary = {
            "status": "complete",
            "excluded_subjects": list(EXCLUDED_SUBJECTS),
            "matching_rule": "exact underscore-delimited subject token in source_sample_id",
            "input_rows": int(matrix.shape[0]),
            "rows_removed": int(remove_mask.sum()),
            "rows_after_removal": int(retained.shape[0]),
            "unique_source_samples_removed": int(excluded["source_sample_id"].nunique()),
            "genes_unchanged": int(matrix.shape[1]),
        }
        selection_summary_path = selection_root / "selection_summary.json"
        safe_write_json(selection_summary_path, selection_summary)

        before_root = branch_root / "02_before_preprocess"
        qc = write_qc_outputs(
            retained,
            retained_metadata,
            before_root,
            title=(
                f"After Ng_MOG1 and Ng_HD3 removal: "
                f"{network.network_id} / {strategy}"
            ),
            include_plotlyjs="cdn",
        )
        preprocess_summary, preprocess_paths = run_diagnostic(
            network,
            strategy,
            retained,
            retained_metadata,
            retained_matrix_path,
        )
        manifest_path = branch_root / "run_manifest.json"
        write_run_manifest(
            manifest_path,
            pipeline_name="transcriptomics_confirmed_outlier_removal_cut_diagnostic",
            run_root=config.run_root,
            parameters={
                "threshold": network.threshold,
                "strategy": strategy,
                "mode": network.mode,
                "specimen": network.specimen,
                "cell_type": network.cell_type,
                "excluded_subjects": list(EXCLUDED_SUBJECTS),
                "cut": "inf",
            },
            generation_rules={
                "input": "existing variance-filtered selected-expression matrix",
                "sample_change_only": True,
                "expression_transform": "none",
                "imputation": "none",
                "preprocess_input_feature_set": (
                    "unchanged from Before Outlier Removal; PyWGCNA preprocess may "
                    "subsequently remove features through its default goodSamplesGenes checks"
                ),
                "preprocess": "PyWGCNA preprocess with cut=inf only",
                "cut_suggestion": "midpoint of largest consecutive sample-linkage height gap",
                "cut_application": "not applied; pending manual review",
            },
            input_paths=[selected_matrix_path, config.clinical_path],
            output_paths=[
                excluded_path,
                retained_path,
                retained_matrix_path,
                selection_summary_path,
                *qc.paths,
                *preprocess_paths,
            ],
            code_paths=[
                Path(__file__),
                Path(__file__).parent / "transcriptomics_qc.py",
                Path(__file__).parent / "run_transcriptomics_variance_filtered.py",
            ],
            packages=(
                "numpy",
                "pandas",
                "matplotlib",
                "scipy",
                "scikit-learn",
                "plotly",
                "PyWGCNA",
            ),
        )
        summary_rows.append(
            {
                "threshold": network.threshold,
                "strategy": strategy,
                "mode": network.mode,
                "specimen": network.specimen,
                "cell_type": network.cell_type,
                "network_id": network.network_id,
                "status": preprocess_summary.get("status"),
                "input_rows_before_removal": int(matrix.shape[0]),
                "rows_removed": int(remove_mask.sum()),
                "rows_after_removal": int(retained.shape[0]),
                "largest_gap_suggested_cut": preprocess_summary.get(
                    "largest_gap_suggested_cut"
                ),
                "downstream_input_status": preprocess_summary.get(
                    "downstream_input_status"
                ),
            }
        )
        del matrix, metadata, retained, retained_metadata
        gc.collect()

    stage_root = (
        config.experiment_root
        / "04_variance_filtered/transcriptomics/outlier_removal_summary"
    )
    stage_root.mkdir(parents=True, exist_ok=True)
    summary_path = stage_root / "outlier_removal_cut_diagnostic_summary.csv"
    new_summary = pd.DataFrame(summary_rows)
    if resume and summary_path.exists():
        old = pd.read_csv(summary_path)
        new_summary = pd.concat([old, new_summary], ignore_index=True).drop_duplicates(
            ["threshold", "strategy", "mode", "specimen", "cell_type"],
            keep="last",
        )
    new_summary.to_csv(summary_path, index=False)
    safe_write_json(
        stage_root / "review_gate.json",
        {
            "status": "proceed_without_finite_sample_cut",
            "excluded_subjects": list(EXCLUDED_SUBJECTS),
            "networks_recorded": int(len(new_summary)),
            "selected_cut_stage": "absent_by_design",
            "next_step": "run findModules from the cut=inf preprocess matrix",
        },
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
            args.threshold
            or ("threshold_5_cells", "threshold_10_cells")
        ),
        only_network=args.only_network,
        resume=not args.no_resume,
    )
    print("Completed confirmed-outlier removal cut=inf preprocess diagnostics.")


if __name__ == "__main__":
    main()
