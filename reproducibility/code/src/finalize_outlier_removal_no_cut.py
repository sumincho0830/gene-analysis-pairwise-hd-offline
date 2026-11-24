"""Finalize the Ng_MOG1/Ng_HD3 branch as cut=inf with no selected-cut stage."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
from typing import Any

import pandas as pd

from run_transcriptomics_stage1 import PipelineConfig
from run_transcriptomics_variance_filtered import STRATEGIES, load_networks


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def safe_remove_tree(path: Path, expected_parent: Path) -> bool:
    if not path.exists():
        return False
    resolved = path.resolve()
    parent = expected_parent.resolve()
    if resolved.parent != parent:
        raise ValueError(f"Refusing to remove unexpected path: {resolved}")
    shutil.rmtree(resolved)
    return True


def run(run_root: Path, offline_root: Path | None) -> dict[str, int]:
    config = PipelineConfig(run_root=run_root)
    networks = load_networks(config, STRATEGIES)
    removed_selected = 0
    summaries_updated = 0
    manifests_updated = 0
    for network in networks:
        branch = network.output_root / "outlier_removal"
        preprocess_parent = branch / "03_preprocess"
        selected = preprocess_parent / "02_selected_cut"
        removed_selected += int(safe_remove_tree(selected, preprocess_parent))

        summary_path = preprocess_parent / "01_cut_diagnostic/preprocess_summary.json"
        if summary_path.exists():
            summary = read_json(summary_path)
            if summary.get("status") == "complete_cut_diagnostic_pending_manual_review":
                summary["status"] = "complete_preprocess_cut_inf"
            summary.pop("selected_cut_status", None)
            summary["finite_cut_applied"] = False
            summary["downstream_input_status"] = (
                "approved_without_finite_sample_cut"
                if summary.get("status") == "complete_preprocess_cut_inf"
                else "not_applicable_preprocess_incomplete"
            )
            write_json(summary_path, summary)
            summaries_updated += 1

        manifest_path = branch / "run_manifest.json"
        if manifest_path.exists():
            manifest = read_json(manifest_path)
            manifest["outputs"] = [
                item
                for item in manifest.get("outputs", [])
                if "/03_preprocess/02_selected_cut/" not in "/" + str(item.get("path", ""))
            ]
            rules = dict(manifest.get("generation_rules", {}))
            rules["sample_tree_cut"] = "none; cut=inf preprocess output proceeds directly to findModules"
            rules["selected_cut_stage"] = "absent by design"
            manifest["generation_rules"] = rules
            manifest["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(manifest_path, manifest)
            manifests_updated += 1

    summary_root = (
        config.experiment_root
        / "04_variance_filtered/transcriptomics/outlier_removal_summary"
    )
    old_summary = summary_root / "outlier_removal_cut_diagnostic_summary.csv"
    new_summary = summary_root / "outlier_removal_preprocess_summary.csv"
    if old_summary.exists():
        frame = pd.read_csv(old_summary)
        if "status" in frame:
            frame["status"] = frame["status"].replace(
                {"complete_cut_diagnostic_pending_manual_review": "complete_preprocess_cut_inf"}
            )
        if "selected_cut_status" in frame:
            frame = frame.drop(columns=["selected_cut_status"])
        frame["finite_cut_applied"] = False
        frame["downstream_input_status"] = frame["status"].map(
            lambda value: (
                "approved_without_finite_sample_cut"
                if value == "complete_preprocess_cut_inf"
                else "not_applicable_preprocess_incomplete"
            )
        )
        frame.to_csv(new_summary, index=False)
        old_summary.unlink()
    write_json(
        summary_root / "review_gate.json",
        {
            "status": "proceed_without_finite_sample_cut",
            "excluded_subjects": ["Ng_MOG1", "Ng_HD3"],
            "networks_recorded": len(networks),
            "selected_cut_stage": "absent_by_design",
            "next_step": "findModules from cut=inf preprocess matrices",
        },
    )

    review_page = (
        config.experiment_root
        / "04_variance_filtered/transcriptomics/outlier_removal_cut_review"
    )
    removed_review_pages = int(
        safe_remove_tree(
            review_page,
            config.experiment_root / "04_variance_filtered/transcriptomics",
        )
    )
    if offline_root is not None:
        offline_page = offline_root / "outlier-removal-cut-review"
        removed_review_pages += int(safe_remove_tree(offline_page, offline_root))
    return {
        "networks": len(networks),
        "selected_cut_directories_removed": removed_selected,
        "preprocess_summaries_updated": summaries_updated,
        "manifests_updated": manifests_updated,
        "review_pages_removed": removed_review_pages,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--offline-root", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(
        json.dumps(
            run(
                args.run_root.resolve(),
                args.offline_root.resolve() if args.offline_root else None,
            ),
            indent=2,
        )
    )
