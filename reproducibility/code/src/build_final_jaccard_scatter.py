"""Package the saved CSF permutation scatter datasets."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import gzip
import io
import json
from pathlib import Path
import shutil


ASSET_DIR = Path("figures/final_jaccard_scatter")
SOURCE_DIR = Path(
    "experiments/05_cross_omics_module_similarity/final_jaccard_before_outlier_plots"
)
BEFORE_CONDITIONS = (
    ("combat_v2_cell_diagnosis_tissue", "ComBat v2", 5, 3600, "combat_v2_5_cells"),
    ("combat_v1_cell_diagnosis", "ComBat v1", 10, 2088, "combat_v1_10_cells"),
)
AFTER_STATISTICS = Path(
    "figures/final_jaccard_bubble/"
    "CSF_threshold_10_cells_hvg_top9000_after_outlier_kME_weighted_Jaccard_permutation_bubble_statistics.csv.gz"
)


@dataclass(frozen=True)
class ScatterSource:
    correction: str
    label: str
    cells: int
    expected: int
    stem: str
    path: Path
    after_outlier: bool = False


def build_final_jaccard_scatter(run_root: Path, dashboard_root: Path) -> int:
    """Place CSV downloads, compact plot coordinates, and the viewer together."""
    before = [
        ScatterSource(correction, label, cells, expected, stem, run_root / SOURCE_DIR /
                      f"{stem}_hvg9000_before_outlier_jaccard_permutation_p.csv")
        for correction, label, cells, expected, stem in BEFORE_CONDITIONS
    ]
    sources = before if any(source.path.is_file() for source in before) else []
    after_path = dashboard_root / AFTER_STATISTICS
    if after_path.is_file():
        sources.append(ScatterSource("combat_v1_cell_diagnosis", "ComBat v1", 10, 3060,
                                     "combat_v1_10_cells", after_path, after_outlier=True))
    if not sources:
        return 0

    index_text = (dashboard_root / "cross-omics-overlap/overlap-index.js").read_text(encoding="utf-8")
    overlap = json.loads(index_text.split("=", 1)[1].strip().removesuffix(";"))["records"]
    datasets = {}
    downloads = {}
    for source in sources:
        csv_bytes = source.path.read_bytes()
        if source.after_outlier:
            csv_bytes = gzip.decompress(csv_bytes)
        filename = source.path.name.removesuffix(".gz")
        downloads[filename] = csv_bytes
        rows = list(csv.DictReader(io.StringIO(csv_bytes.decode("utf-8-sig"), newline="")))
        if len(rows) != source.expected:
            raise ValueError(f"{filename}: expected {source.expected} module pairs, found {len(rows)}")

        records = {}
        saved_pairs = {}
        points = []
        seen = set()
        for row in rows:
            if source.after_outlier:
                record_id = f"{row['record_id']}::after_outlier"
                pair = row["transcript_module"], row["proteomics_module"]
                p_value = float(row["permutation_p_value"])
                q_value = float(row["bh_fdr_q_within_figure"])
                shared = int(row["shared_genes"])
            else:
                record_id = row["rna_record_id"]
                pair = row["rna_module"], row["protein_module"]
                p_value = float(row["permutation_p"])
                q_value = float(row["bh_fdr_q_within_condition"])
                shared = int(row["shared_gene_count"])
            if record_id not in saved_pairs:
                entry = overlap[record_id]
                expected_branch = "after_outlier_removal" if source.after_outlier else "baseline"
                expected = {
                    "analysis_branch": expected_branch, "batch_correction": source.correction,
                    "cell_type": row["cell_type"], "threshold": f"threshold_{source.cells}_cells",
                    "strategy": "02_hvg_top9000", "mode": "cell_separate",
                    "specimen": "CSF", "protein_specimen": "CSF",
                }
                if any(entry.get(field) != value for field, value in expected.items()):
                    raise ValueError(f"Unexpected analysis condition for {record_id}")
                names_for_pairs = entry["pair_columns"]
                rna, protein, score = [names_for_pairs.index(field) for field in (
                    "transcript_module", "proteomics_module", "weighted_jaccard",
                )]
                saved_pairs[record_id] = {(p[rna], p[protein]): p[score] for p in entry["pairs"]}
                records[record_id] = [len(records), row["cell_type"]]
            key = record_id, *pair
            weighted_jaccard = float(row["kme_weighted_jaccard"])
            if (key in seen or int(row["permutations"]) != 10_000
                    or abs(saved_pairs[record_id][pair] - weighted_jaccard) > 1e-10
                    or not 0 < p_value <= q_value <= 1):
                raise ValueError(f"Unexpected scatter row for {key}")
            seen.add(key)
            points.append([records[record_id][0], *pair, weighted_jaccard, p_value, q_value, shared])
        key = f"{source.correction}_after_outlier" if source.after_outlier else source.correction
        datasets[key] = {
            "correction": source.correction, "afterOutlier": source.after_outlier,
            "label": source.label, "cells": source.cells, "stem": source.stem,
            "csv": (ASSET_DIR / filename).as_posix(),
            "records": [[record_id, cell_type] for record_id, (_, cell_type) in records.items()],
            "points": points,
        }

    output = dashboard_root / ASSET_DIR
    output.mkdir(parents=True, exist_ok=True)
    for filename, csv_bytes in downloads.items():
        (output / filename).write_bytes(csv_bytes)
    shutil.copyfile(Path(__file__).with_name("final_jaccard_scatter.js"), output / "viewer-scatter.js")
    (output / "viewer-scatter-data.js").write_text(
        "window.FINAL_JACCARD_SCATTER_DATA="
        + json.dumps(datasets, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        + ";\n",
        encoding="utf-8",
    )
    return len(datasets)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--dashboard-root", type=Path, required=True)
    args = parser.parse_args()
    print(f"Packaged {build_final_jaccard_scatter(args.run_root, args.dashboard_root)} scatter conditions.")
