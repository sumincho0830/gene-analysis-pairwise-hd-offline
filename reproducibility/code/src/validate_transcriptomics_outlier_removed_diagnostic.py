"""Validate the confirmed-outlier-removal diagnostic branch."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

from run_transcriptomics_outlier_removed_diagnostic import (
    EXCLUDED_SUBJECTS,
    excluded_subject,
)
from run_transcriptomics_stage1 import PipelineConfig
from run_transcriptomics_variance_filtered import STRATEGIES, load_networks


QC_FILES = (
    "01_descriptive_statistics/data_summary.csv",
    "01_descriptive_statistics/intensity_distribution.png",
    "01_descriptive_statistics/sample_statistics.csv",
    "01_descriptive_statistics/gene_statistics.csv",
    "01_descriptive_statistics/outlier_candidates.csv",
    "02_data_distribution/heatmap.png",
    "02_data_distribution/pca_metadata.png",
    "02_data_distribution/pca_3d.png",
    "02_data_distribution/pca_3d_interactive.html",
    "02_data_distribution/pca_scores.csv",
)


def validate(run_root: Path) -> dict[str, int]:
    config = PipelineConfig(run_root=run_root)
    networks = load_networks(config, STRATEGIES)
    if len(networks) != 168:
        raise AssertionError(f"Expected 168 networks, found {len(networks)}")
    removed_rows = 0
    networks_with_removed_rows = 0
    suggested_cuts = 0
    cautious_suggestions = 0
    stage_summary_path = (
        config.experiment_root
        / "04_variance_filtered/transcriptomics/outlier_removal_summary"
        / "outlier_removal_preprocess_summary.csv"
    )
    stage_summary = pd.read_csv(stage_summary_path)
    for network in networks:
        branch = network.output_root / "outlier_removal"
        selection = json.loads(
            (branch / "01_sample_selection/selection_summary.json").read_text(
                encoding="utf-8"
            )
        )
        excluded = pd.read_csv(branch / "01_sample_selection/excluded_samples.csv")
        retained = pd.read_csv(branch / "01_sample_selection/retained_samples.csv")
        if set(excluded.get("excluded_subject", [])) - set(EXCLUDED_SUBJECTS):
            raise AssertionError(f"Unexpected exclusion in {network.output_root}")
        if any(
            excluded_subject(value) is not None
            for value in retained["source_sample_id"].astype(str)
        ):
            raise AssertionError(f"Confirmed outlier retained in {network.output_root}")
        if int(selection["rows_removed"]) != len(excluded):
            raise AssertionError(f"Selection count mismatch in {network.output_root}")
        if int(selection["input_rows"]) - int(selection["rows_removed"]) != int(
            selection["rows_after_removal"]
        ):
            raise AssertionError(f"Row arithmetic mismatch in {network.output_root}")
        removed_rows += int(selection["rows_removed"])
        networks_with_removed_rows += int(selection["rows_removed"] > 0)

        qc_root = branch / "02_before_preprocess"
        missing = [relative for relative in QC_FILES if not (qc_root / relative).is_file()]
        if missing:
            raise FileNotFoundError(f"Missing QC outputs in {qc_root}: {missing}")
        pca_scores = pd.read_csv(qc_root / "02_data_distribution/pca_scores.csv")
        id_column = "source_sample_id" if "source_sample_id" in pca_scores else "sample_id"
        if any(
            excluded_subject(value) is not None
            for value in pca_scores[id_column].astype(str)
        ):
            raise AssertionError(f"Confirmed outlier present in PCA: {qc_root}")

        diagnostic = branch / "03_preprocess/01_cut_diagnostic"
        summary = json.loads(
            (diagnostic / "preprocess_summary.json").read_text(encoding="utf-8")
        )
        recommendation_path = diagnostic / "cut_recommendation.json"
        recommendation = (
            json.loads(recommendation_path.read_text(encoding="utf-8"))
            if recommendation_path.exists()
            else {"applied": False, "suggested_cut": None}
        )
        finite_cut_applied = summary.get(
            "finite_cut_applied", summary.get("cut_applied")
        )
        if finite_cut_applied is not False:
            raise AssertionError(f"Finite cut applied in {diagnostic}")
        if recommendation.get("applied") is not False:
            raise AssertionError(f"Recommendation applied in {diagnostic}")
        suggested_cut = recommendation.get("suggested_cut")
        suggested_cuts += int(suggested_cut is not None)
        cautious_suggestions += int(bool(recommendation.get("requires_caution")))
        if suggested_cut is not None:
            candidates = pd.read_csv(diagnostic / "cut_candidates.csv")
            selected_rows = candidates.loc[
                candidates["selected_review_suggestion"].astype(bool)
            ]
            if len(selected_rows) != 1:
                raise AssertionError(f"Expected one selected gap in {diagnostic}")
            selected_row = selected_rows.iloc[0]
            maximum_gap = float(candidates["absolute_gap"].max())
            if not math.isclose(
                float(selected_row["absolute_gap"]), maximum_gap, rel_tol=1e-12
            ):
                raise AssertionError(f"Suggestion is not the global largest gap: {diagnostic}")
            midpoint = (
                float(selected_row["lower_height"])
                + float(selected_row["upper_height"])
            ) / 2.0
            if not math.isclose(float(suggested_cut), midpoint, rel_tol=1e-12):
                raise AssertionError(f"Suggestion is not the gap midpoint: {diagnostic}")

            strategy = next(
                value for value in STRATEGIES if value in network.output_root.parts
            )
            stage_row = stage_summary.loc[
                (stage_summary["threshold"] == network.threshold)
                & (stage_summary["strategy"] == strategy)
                & (stage_summary["mode"] == network.mode)
                & (stage_summary["specimen"] == network.specimen)
                & (stage_summary["cell_type"] == network.cell_type)
            ]
            if len(stage_row) != 1 or not math.isclose(
                float(stage_row.iloc[0]["largest_gap_suggested_cut"]),
                float(suggested_cut),
                rel_tol=1e-12,
            ):
                raise AssertionError(f"Stage summary is stale for {diagnostic}")
        selected = branch / "03_preprocess/02_selected_cut"
        if selected.exists():
            raise AssertionError(f"Selected-cut stage must be absent: {selected}")
        if summary.get("status") not in {
            "complete_preprocess_cut_inf",
            "not_run_insufficient_samples",
        }:
            raise AssertionError(f"Unexpected cut=inf preprocess status in {diagnostic}")

    return {
        "networks": len(networks),
        "networks_with_removed_rows": networks_with_removed_rows,
        "removed_matrix_rows": removed_rows,
        "networks_with_cut_suggestion": suggested_cuts,
        "suggestions_requiring_caution": cautious_suggestions,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(json.dumps(validate(args.run_root.resolve()), indent=2))
