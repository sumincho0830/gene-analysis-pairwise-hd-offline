"""Add confirmed-outlier-removal artifacts to an existing offline release.

This is intentionally an incremental packager: it preserves the current
dashboard, review state schema, module assets, and hand-tuned UI while adding
the new transcriptomics branch generated under each variance-filtered network.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any

from refresh_pca_assets import (
    DATA_PREFIX,
    convert_lossless_webp,
    read_bundle,
    rebuild_checksums,
    sha256,
    write_bundle,
)


ARTIFACT_SPECS = (
    (
        "after_outlier_heatmap",
        "Heatmap",
        "outlier_removal/02_before_preprocess/02_data_distribution/heatmap.png",
        "image",
        "before preprocess",
    ),
    (
        "after_outlier_pca",
        "PCA metadata",
        "outlier_removal/02_before_preprocess/02_data_distribution/pca_metadata.png",
        "image",
        "before preprocess",
    ),
    (
        "after_outlier_pca_3d",
        "PCA 3D interactive",
        "outlier_removal/02_before_preprocess/02_data_distribution/pca_3d_interactive.html",
        "html",
        "before preprocess",
    ),
    (
        "after_outlier_cut_diagnostic",
        "Cut diagnostic",
        "outlier_removal/03_preprocess/01_cut_diagnostic/dendrogram_cut_diagnostic.png",
        "image",
        "preprocess",
    ),
    (
        "after_outlier_after_heatmap",
        "After heatmap",
        "outlier_removal/04_after_preprocess/02_data_distribution/heatmap.png",
        "image",
        "after preprocess after outlier removal",
    ),
    (
        "after_outlier_after_pca",
        "After PCA",
        "outlier_removal/04_after_preprocess/02_data_distribution/pca_metadata.png",
        "image",
        "after preprocess after outlier removal",
    ),
    (
        "after_outlier_after_pca_3d",
        "After PCA 3D interactive",
        "outlier_removal/04_after_preprocess/02_data_distribution/pca_3d_interactive.html",
        "html",
        "after preprocess after outlier removal",
    ),
    (
        "after_outlier_selected_cut",
        "Selected cut",
        "outlier_removal/03_preprocess/02_selected_cut/selected_cut_pending.html",
        "html",
        "preprocess",
    ),
)

FILE_SPECS = (
    (
        "outlier_removal_excluded_samples",
        "Confirmed excluded samples",
        "outlier_removal/01_sample_selection/excluded_samples.csv",
        "outlier removal",
    ),
    (
        "outlier_removal_selection_summary",
        "Outlier-removal selection summary",
        "outlier_removal/01_sample_selection/selection_summary.json",
        "outlier removal",
    ),
    (
        "after_outlier_pca_scores",
        "After-outlier-removal PCA scores",
        "outlier_removal/02_before_preprocess/02_data_distribution/pca_scores.csv",
        "outlier removal QC",
    ),
    (
        "after_outlier_after_pca_scores",
        "After-outlier post-preprocess PCA scores",
        "outlier_removal/04_after_preprocess/02_data_distribution/pca_scores.csv",
        "after preprocess after outlier removal",
    ),
    (
        "after_outlier_cut_recommendation",
        "Largest-gap cut recommendation",
        "outlier_removal/03_preprocess/01_cut_diagnostic/cut_recommendation.json",
        "outlier removal preprocess",
    ),
    (
        "after_outlier_cut_candidates",
        "All cut candidates",
        "outlier_removal/03_preprocess/01_cut_diagnostic/cut_candidates.csv",
        "outlier removal preprocess",
    ),
    (
        "after_outlier_preprocess_summary",
        "Cut-diagnostic preprocess summary",
        "outlier_removal/03_preprocess/01_cut_diagnostic/preprocess_summary.json",
        "outlier removal preprocess",
    ),
    (
        "after_outlier_selected_cut_status",
        "Selected-cut review status",
        "outlier_removal/03_preprocess/02_selected_cut/review_status.json",
        "outlier removal preprocess",
    ),
    (
        "after_outlier_run_manifest",
        "Outlier-removal run manifest",
        "outlier_removal/run_manifest.json",
        "provenance",
    ),
)


def network_root(source_run: Path, record: dict[str, Any]) -> Path:
    dimensions = record["dimensions"]
    root = (
        source_run
        / "experiments/04_variance_filtered/transcriptomics"
        / dimensions["threshold"]
        / dimensions["strategy"]
        / dimensions["mode"]
        / dimensions["specimen"]
    )
    if dimensions["mode"] == "cell_separate":
        root /= dimensions["cell_type"]
    return root


def replace_plotly_cdn(html: str) -> str:
    replacement = '<script charset="utf-8" src="../vendor/plotly.min.js"></script>'
    return re.sub(
        r'<script[^>]+src="https://cdn\.plot\.ly/[^\"]+"[^>]*></script>',
        replacement,
        html,
    )


def copy_html(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    content = replace_plotly_cdn(source.read_text(encoding="utf-8"))
    target.write_text(content, encoding="utf-8")


def replace_by_id(items: list[dict[str, Any]], item: dict[str, Any]) -> None:
    items[:] = [existing for existing in items if existing.get("id") != item["id"]]
    items.append(item)


def refresh(
    source_run: Path,
    release_root: Path,
    workers: int,
) -> dict[str, Any]:
    bundle_path = release_root / "dashboard-data.js"
    bundle = read_bundle(bundle_path)
    conversions: dict[str, tuple[Path, Path]] = {}
    copied_html = 0
    copied_data = 0
    records_updated = 0
    missing: list[str] = []

    for record in bundle["records"]:
        if record.get("source_id") != "transcriptomics_variance_findmodules":
            continue
        # The removal-state tabs now carry the before/after meaning, so keep
        # the leaf artifact labels identical for direct visual comparison.
        baseline_labels = {
            "before_heatmap": "Heatmap",
            "before_pca": "PCA metadata",
            "before_pca_3d": "PCA 3D interactive",
        }
        for artifact in record.get("artifacts", []):
            if artifact.get("id") in baseline_labels:
                artifact["label"] = baseline_labels[artifact["id"]]
                artifact["outlier_removal_state"] = "before"
        root = network_root(source_run, record)
        record_hash = hashlib.sha256(record["id"].encode("utf-8")).hexdigest()[:16]
        added = 0
        for artifact_id, label, relative, kind, stage in ARTIFACT_SPECS:
            source = root / relative
            if not source.is_file():
                missing.append(f"{record['id']}::{relative}")
                continue
            if kind == "image":
                digest = sha256(source)[:16]
                target = release_root / "assets" / artifact_id / f"{digest}.webp"
                conversions[target.as_posix().lower()] = (source, target)
            else:
                digest = sha256(source)[:16]
                target = release_root / "assets" / artifact_id / f"{digest}.html"
                copy_html(source, target)
                copied_html += 1
            replace_by_id(
                record.setdefault("artifacts", []),
                {
                    "id": artifact_id,
                    "label": label,
                    "stage": stage,
                    "kind": kind,
                    "path": target.relative_to(release_root).as_posix(),
                    "outlier_removal_state": "after",
                },
            )
            added += 1

        for file_id, label, relative, category in FILE_SPECS:
            source = root / relative
            if not source.is_file():
                missing.append(f"{record['id']}::{relative}")
                continue
            target_name = Path(relative).name
            # Both QC stages contain a pca_scores.csv. Preserve the existing
            # before-preprocess path and give the post-preprocess table its own
            # filename so one stage cannot overwrite the other.
            if file_id == "after_outlier_after_pca_scores":
                target_name = "after_outlier_after_pca_scores.csv"
            target = (
                release_root
                / "data/outlier_removal"
                / record_hash
                / target_name
            )
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            copied_data += 1
            replace_by_id(
                record.setdefault("files", []),
                {
                    "id": file_id,
                    "label": label,
                    "category": category,
                    "path": target.relative_to(release_root).as_posix(),
                },
            )
        if added:
            records_updated += 1

    tasks = [
        (str(source), str(target))
        for source, target in conversions.values()
        if not target.is_file()
    ]
    converted = 0
    with ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = [executor.submit(convert_lossless_webp, task) for task in tasks]
        for future in as_completed(futures):
            future.result()
            converted += 1
            if converted % 25 == 0 or converted == len(futures):
                print(
                    f"Converted {converted}/{len(futures)} outlier-removal assets",
                    flush=True,
                )

    now = datetime.now(timezone.utc).isoformat()
    bundle["generated_at_utc"] = now
    write_bundle(bundle_path, bundle)
    manifest_path = release_root / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outlier_removal_refresh"] = {
        "generated_at_utc": now,
        "excluded_subjects": ["Ng_MOG1", "Ng_HD3"],
        "records_updated": records_updated,
        "image_assets": len(conversions),
        "html_assets": copied_html,
        "data_files": copied_data,
        "missing_sources": len(missing),
        "selected_cut_status": "pending_manual_review",
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    manifest["release_files"] = rebuild_checksums(release_root)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    rebuild_checksums(release_root)
    return {
        "records_updated": records_updated,
        "unique_image_assets": len(conversions),
        "converted_image_assets": converted,
        "copied_html": copied_html,
        "copied_data": copied_data,
        "missing": missing,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument(
        "--release-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = refresh(
        args.source_run.resolve(),
        args.release_root.resolve(),
        args.workers,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
