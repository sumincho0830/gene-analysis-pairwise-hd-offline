"""Validate the no-finite-cut transcriptomics downstream branch."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


CONDITION = "01_default_ds2_mms50_med0p20"
EXCLUDED_SUBJECTS = ("Ng_MOG1", "Ng_HD3")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def network_root(stage: Path, row: Any) -> Path:
    root = stage / row.threshold / row.strategy / row.mode / row.specimen
    if row.mode == "cell_separate":
        root /= row.cell_type
    return root


def assert_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing or empty file: {path}")


def run(run_root: Path) -> dict[str, Any]:
    stage = run_root / "experiments/04_variance_filtered/transcriptomics"
    preprocess_path = (
        stage
        / "outlier_removal_summary/outlier_removal_preprocess_summary.csv"
    )
    find_path = stage / "outlier_removal_find_modules_summary.csv"
    downstream_path = stage / "outlier_removal_module_diagnosis_summary.csv"
    cross_root = (
        run_root
        / "experiments/05_cross_omics_module_similarity/after_outlier_removal"
    )
    cross_path = cross_root / "cross_omics_similarity_summary.csv"
    for path in (preprocess_path, find_path, downstream_path, cross_path):
        assert_file(path)

    preprocess = pd.read_csv(preprocess_path)
    find = pd.read_csv(find_path)
    downstream = pd.read_csv(downstream_path)
    cross = pd.read_csv(cross_path)
    if len(preprocess) != 168 or len(find) != 168:
        raise AssertionError(
            f"Expected 168 preprocess/find rows, got {len(preprocess)}/{len(find)}"
        )
    ready = preprocess["status"].eq("complete_preprocess_cut_inf").sum()
    insufficient = preprocess["status"].eq("not_run_insufficient_samples").sum()
    if (ready, insufficient) != (165, 3):
        raise AssertionError(f"Unexpected preprocess status counts: {ready}, {insufficient}")
    complete = find.loc[find["status"].eq("complete_findModules")].copy()
    skipped = find.loc[find["status"].eq("not_run_preprocess_incomplete")].copy()
    failures = find.loc[find["status"].str.contains("fail|error", case=False, na=False)]
    if len(complete) != ready or len(skipped) != insufficient or not failures.empty:
        raise AssertionError(
            f"Unexpected findModules counts: complete={len(complete)}, "
            f"skipped={len(skipped)}, failed={len(failures)}"
        )
    if complete["finite_sample_cut_applied"].astype(str).str.lower().ne("false").any():
        raise AssertionError("A completed network reports a finite sample cut")
    if complete["feature_filtering_before_findModules"].ne("none").any():
        raise AssertionError("Unexpected feature filtering before findModules")
    if complete["expression_transform_before_findModules"].ne("none").any():
        raise AssertionError("Unexpected expression transform before findModules")

    required_find = (
        "analysis_summary.json",
        "module_assignments.csv.gz",
        "module_eigengenes.csv",
        "module_summary.csv",
        "sample_info.csv",
        "figures/soft_threshold.png",
        "figures/gene_dendrogram_modules.png",
        "figures/module_eigengene_dendrogram.png",
        "figures/module_sizes.png",
    )
    required_downstream = (
        "analysis_summary.json",
        "pairwise_hd_analysis_summary.json",
        "feature_module_kme.csv.gz",
        "module_kme_summary.csv",
        "module_diagnosis_statistics.csv",
        "module_diagnosis_pairwise_hd_statistics.csv",
        "figures/module_diagnosis_heatmap.png",
        "figures/module_diagnosis_pairwise_hd_heatmap.png",
    )
    checked_samples = 0
    for row in complete.itertuples(index=False):
        root = network_root(stage, row)
        selected_cut = root / "outlier_removal/03_preprocess/02_selected_cut"
        if selected_cut.exists():
            raise AssertionError(f"Selected-cut stage must be absent: {selected_cut}")
        find_root = root / "outlier_removal/05_find_modules" / CONDITION
        for relative in required_find:
            assert_file(find_root / relative)
        sample_info = pd.read_csv(find_root / "sample_info.csv", index_col=0)
        sample_ids = sample_info.index.astype(str)
        for subject in EXCLUDED_SUBJECTS:
            if sample_ids.str.contains(subject, regex=False).any():
                raise AssertionError(f"Excluded subject {subject} remains in {find_root}")
        checked_samples += len(sample_info)
        diagnosis_root = root / "outlier_removal/06_kme_module_diagnosis" / CONDITION
        for relative in required_downstream:
            assert_file(diagnosis_root / relative)
        official = load_json(diagnosis_root / "analysis_summary.json")
        pairwise = load_json(diagnosis_root / "pairwise_hd_analysis_summary.json")
        if official.get("status") != "complete" or pairwise.get("status") != "complete":
            raise AssertionError(f"Incomplete diagnosis result: {diagnosis_root}")
        if official.get("sample_tree_cut") != "none_cut_inf":
            raise AssertionError(f"Unexpected sample-tree cut: {diagnosis_root}")

    if len(downstream) != ready:
        raise AssertionError(f"Expected {ready} downstream rows, got {len(downstream)}")
    if not downstream["official_status"].eq("complete").all():
        raise AssertionError("At least one official kME/diagnosis result is incomplete")
    if not downstream["pairwise_status"].eq("complete").all():
        raise AssertionError("At least one pairwise-HD result is incomplete")
    if len(cross) != ready:
        raise AssertionError(f"Expected {ready} cross-omics rows, got {len(cross)}")
    if cross["analysis_branch"].ne("after_outlier_removal").any():
        raise AssertionError("Cross-omics summary contains the wrong branch")
    if cross["maximum_weighted_jaccard"].isna().any():
        raise AssertionError("Cross-omics similarity is missing")

    result = {
        "status": "passed",
        "preprocess_networks": len(preprocess),
        "preprocess_complete_cut_inf": int(ready),
        "preprocess_insufficient_samples": int(insufficient),
        "findModules_complete": len(complete),
        "kme_module_diagnosis_complete": len(downstream),
        "cross_omics_complete": len(cross),
        "finite_sample_cuts_applied": 0,
        "excluded_subjects": list(EXCLUDED_SUBJECTS),
        "sample_rows_checked": checked_samples,
    }
    output = stage / "outlier_removal_downstream_validation.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(json.dumps(run(args.run_root.resolve()), indent=2))
