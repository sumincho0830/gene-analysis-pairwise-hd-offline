"""Package PCA outlier evidence for the offline inspection dashboard."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
from typing import Any

from refresh_pca_assets import rebuild_checksums


DATA_PREFIX = "window.OUTLIER_CANDIDATE_DATA = "
CSV_NAMES = (
    "outlier_candidate_occurrences.csv",
    "outlier_candidate_sample_counts.csv",
    "outlier_candidate_experiment_counts.csv",
)
TRUE_VALUES = {"true", "1", "yes"}
SEVERITY_PERCENTILE_BASELINE = 95.0
CELL_TYPE_ORDER = (
    "B_cell",
    "CD4_T_cell",
    "CD8_T_cell",
    "NK_cell",
    "MAIT_gDT_cell",
    "Classical_Monocytes",
    "Non_classical_Monocytes",
    "Macrophages",
    "cDCs",
    "pDCs",
    "Plasmablast",
    "Proliferating_T_cell",
    "other",
)


def number(value: str) -> float | int | None:
    if value is None or value == "":
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    if not math.isfinite(result):
        return None
    return int(result) if result.is_integer() else result


def dashboard_target(path: str) -> tuple[str | None, str | None]:
    parts = PurePosixPath(path).parts
    if len(parts) < 7 or parts[0] != "04_variance_filtered":
        return None, None
    threshold, strategy, mode, specimen = parts[2:6]
    if mode == "cell_separate":
        cell_type = parts[6]
        tail = parts[7:]
    elif mode == "cell_integrated":
        cell_type = "cell_integrated"
        tail = parts[6:]
    else:
        return None, None
    phase = next(
        (part for part in tail if part in {"02_before_preprocess", "04_after_preprocess"}),
        None,
    )
    artifact_id = {
        "02_before_preprocess": "before_pca",
        "04_after_preprocess": "after_pca",
    }.get(phase)
    record_id = (
        "transcriptomics_variance::"
        f"{threshold}::{strategy}::{mode}::{specimen}::{cell_type}"
    )
    return record_id, artifact_id


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_dimensions(relative: Path) -> dict[str, str]:
    parts = relative.parts[:-1]
    threshold = next((part for part in parts if part.startswith("threshold_")), "")
    mode = next((part for part in parts if part in {"cell_separate", "cell_integrated"}), "")
    specimen = next((part for part in reversed(parts) if part in {"CSF", "PB", "ALL"}), "")
    strategy = next(
        (part for part in parts if re.match(r"^0[1-3]_(?:cpm|hvg)", part)),
        "",
    )
    phase = next(
        (
            part
            for part in parts
            if part
            in {
                "01_descriptive_statistics",
                "02_data_distribution",
                "02_before_preprocess",
                "03_preprocess",
                "04_after_preprocess",
            }
        ),
        "",
    )
    return {
        "stage": parts[0] if parts else "",
        "modality": parts[1] if len(parts) > 1 else "",
        "threshold": threshold,
        "mode": mode,
        "strategy": strategy,
        "specimen": specimen or "ALL",
        "phase": phase,
    }


def cell_sort_key(value: str) -> tuple[int, str]:
    try:
        return CELL_TYPE_ORDER.index(value), value
    except ValueError:
        return len(CELL_TYPE_ORDER), value


def severity_weight(percentile: float | int | None) -> float:
    """Map the upper 5% percentile tail linearly to a bounded 0-1 weight."""

    if not isinstance(percentile, (int, float)):
        return 0.0
    return max(
        0.0,
        min(
            1.0,
            (float(percentile) - SEVERITY_PERCENTILE_BASELINE)
            / (100.0 - SEVERITY_PERCENTILE_BASELINE),
        ),
    )


def build_cell_type_experiment_sets(experiments_root: Path) -> list[dict[str, Any]]:
    """Pivot saved cell-separated PCA results into candidate-by-cell matrices."""

    groups: dict[tuple[str, ...], dict[str, Any]] = {}
    for score_path in sorted(experiments_root.rglob("pca_scores.csv")):
        relative = score_path.relative_to(experiments_root)
        dimensions = parse_dimensions(relative)
        if dimensions["modality"] != "transcriptomics" or dimensions["mode"] != "cell_separate":
            continue
        score_rows = read_csv(score_path)
        if not score_rows:
            continue
        columns = set(score_rows[0])
        distance_metric = (
            "pca_robust_distance_pc1_pc5"
            if "pca_robust_distance_pc1_pc5" in columns
            else "pca_distance_robust_z"
        )
        pca_candidate_column = (
            "pca_outlier_candidate"
            if "pca_outlier_candidate" in columns
            else "pca_outlier"
            if "pca_outlier" in columns
            else ""
        )
        heatmap_candidate_column = (
            "outlier_candidate" if "outlier_candidate" in columns else pca_candidate_column
        )
        distances = [number(row.get(distance_metric, "")) for row in score_rows]
        valid_distances = [float(value) for value in distances if isinstance(value, (int, float))]
        if not valid_distances:
            continue
        group_key = tuple(
            dimensions[name]
            for name in ("stage", "threshold", "strategy", "specimen", "phase", "mode")
        )
        group = groups.setdefault(
            group_key,
            {
                "dimensions": dimensions,
                "cell_types": set(),
                "samples": {},
                "score_files": set(),
            },
        )
        group["score_files"].add(relative.as_posix())
        record_id, artifact_id = dashboard_target(relative.as_posix())
        for row, distance in zip(score_rows, distances):
            if not isinstance(distance, (int, float)):
                continue
            cell_type = row.get("cell_type", "")
            if not cell_type:
                continue
            group["cell_types"].add(cell_type)
            sample_id = row.get("source_sample_id") or row.get("sample_id", "")
            if not sample_id:
                continue
            percentile = 100.0 * sum(value <= float(distance) for value in valid_distances) / len(valid_distances)
            pca_flag = (
                1
                if pca_candidate_column
                and row.get(pca_candidate_column, "").strip().lower() in TRUE_VALUES
                else 0
            )
            heatmap_flag = (
                1
                if heatmap_candidate_column
                and row.get(heatmap_candidate_column, "").strip().lower() in TRUE_VALUES
                else pca_flag
            )
            sample = group["samples"].setdefault(
                sample_id,
                {
                    "source_sample_id": sample_id,
                    "specimen": row.get("specimen", ""),
                    "cohort": row.get("cohort", ""),
                    "diagnosis": row.get("diagnosis", ""),
                    "cells": {},
                },
            )
            cell_result = {
                "flag": pca_flag,
                "pca_flag": pca_flag,
                "heatmap_flag": heatmap_flag,
                "candidate_source": (
                    "pca" if pca_flag else "descriptive" if heatmap_flag else "none"
                ),
                "distance": distance,
                "percentile": percentile,
                "distance_metric": distance_metric,
                "score_path": relative.as_posix(),
                "record_id": record_id,
                "artifact_id": artifact_id,
            }
            previous = sample["cells"].get(cell_type)
            if (
                previous is None
                or pca_flag > previous["pca_flag"]
                or heatmap_flag > previous["heatmap_flag"]
                or percentile > previous["percentile"]
            ):
                sample["cells"][cell_type] = cell_result

    output: list[dict[str, Any]] = []
    for group_key, group in groups.items():
        cell_types = sorted(group["cell_types"], key=cell_sort_key)
        pca_samples = []
        heatmap_samples = []
        for sample in group["samples"].values():
            pca_flagged = sum(cell["pca_flag"] for cell in sample["cells"].values())
            heatmap_flagged = sum(cell["heatmap_flag"] for cell in sample["cells"].values())
            tested = len(sample["cells"])
            pca_weighted = sum(
                severity_weight(cell["percentile"])
                for cell in sample["cells"].values()
                if cell["pca_flag"]
            )
            heatmap_weighted = sum(
                severity_weight(cell["percentile"])
                for cell in sample["cells"].values()
                if cell["heatmap_flag"]
            )
            enriched = {
                **sample,
                "outlier_cell_count": pca_flagged,
                "tested_cell_count": tested,
                "outlier_cell_fraction": pca_flagged / tested if tested else None,
                "pca_outlier_cell_count": pca_flagged,
                "pca_weighted_count": pca_weighted,
                "heatmap_outlier_cell_count": heatmap_flagged,
                "heatmap_weighted_count": heatmap_weighted,
            }
            if pca_flagged:
                pca_samples.append(enriched)
            if heatmap_flagged:
                heatmap_samples.append(enriched)
        if not pca_samples and not heatmap_samples:
            continue
        pca_samples.sort(
            key=lambda row: (
                -row["pca_weighted_count"],
                -row["outlier_cell_count"],
                -row["outlier_cell_fraction"],
                row["source_sample_id"],
            )
        )
        heatmap_samples.sort(
            key=lambda row: (
                -row["heatmap_weighted_count"],
                -row["heatmap_outlier_cell_count"],
                row["source_sample_id"],
            )
        )
        dimensions = group["dimensions"]
        label_parts = [
            dimensions["stage"],
            dimensions["threshold"],
            dimensions["strategy"],
            dimensions["specimen"],
            dimensions["phase"],
        ]
        output.append(
            {
                "id": "::".join(group_key),
                "label": " · ".join(part.replace("_", " ") for part in label_parts if part),
                "dimensions": dimensions,
                "cell_types": cell_types,
                "candidate_samples": pca_samples,
                "heatmap_candidate_samples": heatmap_samples,
                "score_file_count": len(group["score_files"]),
            }
        )
    output.sort(
        key=lambda group: (
            group["dimensions"]["stage"],
            group["dimensions"]["threshold"],
            group["dimensions"]["strategy"],
            group["dimensions"]["specimen"],
            group["dimensions"]["phase"],
        )
    )
    return output


def clean_occurrence(
    row: dict[str, str],
    *,
    distance_percentile: float | None,
    distance_sample_count: int,
) -> dict[str, Any]:
    record_id, artifact_id = dashboard_target(row["experiment_relative_path"])
    distance = number(row.get("pca_distance", ""))
    return {
        "experiment_path": row["experiment_relative_path"],
        "stage": row.get("experiment_stage", ""),
        "modality": row.get("modality", ""),
        "experiment_threshold": row.get("threshold", ""),
        "mode": row.get("mode", ""),
        "strategy": row.get("strategy", ""),
        "path_specimen": row.get("path_specimen", ""),
        "phase": row.get("phase", ""),
        "sample_id": row.get("sample_id", ""),
        "source_sample_id": row.get("source_sample_id", ""),
        "cohort": row.get("cohort", ""),
        "specimen": row.get("specimen", ""),
        "cell_type": row.get("cell_type", ""),
        "diagnosis": row.get("diagnosis", ""),
        "sex": row.get("sex", ""),
        "age": number(row.get("age", "")),
        "candidate_source": row.get("candidate_source", ""),
        "distance_metric": row.get("pca_distance_metric", ""),
        "distance": distance,
        "distance_percentile": distance_percentile,
        "distance_sample_count": distance_sample_count,
        "pc1": number(row.get("PC1", "")),
        "pc2": number(row.get("PC2", "")),
        "pc3": number(row.get("PC3", "")),
        "record_id": record_id,
        "artifact_id": artifact_id,
    }


def add_empirical_percentiles(
    rows: list[dict[str, str]],
    experiments_root: Path,
) -> list[tuple[float | None, int]]:
    """Return ECDF percentile within each complete experiment distance vector."""

    cache: dict[tuple[str, str], list[float]] = {}
    results: list[tuple[float | None, int]] = []
    for row in rows:
        relative = row["experiment_relative_path"]
        metric = row.get("pca_distance_metric", "")
        key = (relative, metric)
        if key not in cache:
            score_path = experiments_root / PurePosixPath(relative)
            scores = []
            if score_path.is_file() and metric:
                for score_row in read_csv(score_path):
                    value = number(score_row.get(metric, ""))
                    if isinstance(value, (int, float)):
                        scores.append(float(value))
            cache[key] = scores
        scores = cache[key]
        distance = number(row.get("pca_distance", ""))
        if not scores or not isinstance(distance, (int, float)):
            results.append((None, len(scores)))
            continue
        percentile = 100.0 * sum(value <= float(distance) for value in scores) / len(scores)
        results.append((percentile, len(scores)))
    return results


def clean_sample(row: dict[str, str]) -> dict[str, Any]:
    return {
        "occurrence_rank": number(row.get("occurrence_rank", "")),
        "source_sample_id": row["source_sample_id"],
        "occurrence_count": number(row.get("occurrence_count", "")),
        "distinct_experiment_count": number(row.get("distinct_experiment_count", "")),
        "distinct_analysis_row_count": number(row.get("distinct_analysis_row_count", "")),
        "modalities": row.get("modalities", ""),
        "specimens": row.get("specimens", ""),
        "cohorts": row.get("cohorts", ""),
        "cell_types": row.get("cell_types", ""),
        "diagnoses": row.get("diagnoses", ""),
        "candidate_sources": row.get("candidate_sources", ""),
    }


def clean_experiment(row: dict[str, str]) -> dict[str, Any]:
    return {
        "experiment_path": row["experiment_relative_path"],
        "candidate_count": number(row.get("candidate_count", "")),
        "distinct_source_sample_count": number(row.get("distinct_source_sample_count", "")),
        "maximum_pca_distance": number(row.get("maximum_pca_distance", "")),
        "distance_metric": row.get("pca_distance_metric", ""),
    }


def build(source_dir: Path, release_root: Path) -> dict[str, Any]:
    raw_occurrences = read_csv(source_dir / "outlier_candidate_occurrences.csv")
    percentile_results = add_empirical_percentiles(raw_occurrences, source_dir.parent)
    occurrences = [
        clean_occurrence(
            row,
            distance_percentile=percentile,
            distance_sample_count=sample_count,
        )
        for row, (percentile, sample_count) in zip(raw_occurrences, percentile_results)
    ]
    samples = [
        clean_sample(row)
        for row in read_csv(source_dir / "outlier_candidate_sample_counts.csv")
    ]
    experiments = [
        clean_experiment(row)
        for row in read_csv(source_dir / "outlier_candidate_experiment_counts.csv")
    ]
    cell_type_experiment_sets = build_cell_type_experiment_sets(source_dir.parent)
    summary = json.loads((source_dir / "summary.json").read_text(encoding="utf-8"))
    output_dir = release_root / "outlier-candidates"
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = release_root / "data" / "QC" / "outlier_candidates"
    data_dir.mkdir(parents=True, exist_ok=True)
    for name in CSV_NAMES[1:]:
        shutil.copy2(source_dir / name, data_dir / name)
    enriched_occurrence_path = data_dir / CSV_NAMES[0]
    occurrence_fields = list(raw_occurrences[0]) + [
        "pca_distance_percentile",
        "pca_distance_sample_count",
    ]
    with enriched_occurrence_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=occurrence_fields)
        writer.writeheader()
        for row, (percentile, sample_count) in zip(raw_occurrences, percentile_results):
            writer.writerow(
                {
                    **row,
                    "pca_distance_percentile": percentile,
                    "pca_distance_sample_count": sample_count,
                }
            )
    shutil.copy2(source_dir / "summary.json", data_dir / "summary.json")

    bundle = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_rule": summary.get("candidate_rule"),
        "severity_weight": {
            "name": "upper-tail percentile weight",
            "formula": "max(0, min(1, (percentile - 95) / 5))",
            "percentile_baseline": SEVERITY_PERCENTILE_BASELINE,
            "purpose": "review prioritization only; raw counts and saved candidate flags remain unchanged",
        },
        "summary": summary,
        "samples": samples,
        "experiments": experiments,
        "occurrences": occurrences,
        "cell_type_experiment_sets": cell_type_experiment_sets,
        "downloads": {
            "occurrences": "data/QC/outlier_candidates/outlier_candidate_occurrences.csv",
            "samples": "data/QC/outlier_candidates/outlier_candidate_sample_counts.csv",
            "experiments": "data/QC/outlier_candidates/outlier_candidate_experiment_counts.csv",
        },
    }
    (output_dir / "outlier-candidates.js").write_text(
        DATA_PREFIX + json.dumps(bundle, ensure_ascii=False, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )

    manifest_path = release_root / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outlier_review"] = {
        "candidate_occurrences": len(occurrences),
        "candidate_samples": len(samples),
        "experiments_with_candidates": len(experiments),
        "generated_at_utc": bundle["generated_at_utc"],
        "distance_percentile_method": "empirical CDF within each complete experiment distance vector",
        "cell_type_experiment_sets": len(cell_type_experiment_sets),
        "severity_weight_formula": "max(0, min(1, (percentile - 95) / 5))",
        "severity_weight_role": "review prioritization only; does not replace raw candidate counts",
    }
    manifest["release_files"] = sum(
        1
        for path in release_root.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS.txt"
        and ".git" not in path.relative_to(release_root).parts
        and "__pycache__" not in path.relative_to(release_root).parts
        and path.suffix.lower() != ".pyc"
    ) + 1
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    rebuild_checksums(release_root)
    return manifest["outlier_review"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument(
        "--release-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build(args.source_dir.resolve(), args.release_root.resolve())
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
