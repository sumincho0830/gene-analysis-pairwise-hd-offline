"""Package the two saved before-outlier permutation scatter datasets."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil


ASSET_DIR = Path("figures/final_jaccard_scatter")
SOURCE_DIR = Path(
    "experiments/05_cross_omics_module_similarity/final_jaccard_before_outlier_plots"
)
CONDITIONS = (
    ("combat_v2_cell_diagnosis_tissue", "ComBat v2", 5, 3600, "combat_v2_5_cells"),
    ("combat_v1_cell_diagnosis", "ComBat v1", 10, 2088, "combat_v1_10_cells"),
)


def build_final_jaccard_scatter(run_root: Path, dashboard_root: Path) -> int:
    """Place CSV downloads, compact plot coordinates, and the viewer together."""
    names = [f"{stem}_hvg9000_before_outlier_jaccard_permutation_p.csv"
             for _, _, _, _, stem in CONDITIONS]
    sources = [run_root / SOURCE_DIR / name for name in names]
    if not any(source.is_file() for source in sources):
        return 0

    index_text = (dashboard_root / "cross-omics-overlap/overlap-index.js").read_text(encoding="utf-8")
    overlap = json.loads(index_text.split("=", 1)[1].strip().removesuffix(";"))["records"]
    datasets = {}
    for (correction, label, cells, expected, stem), source in zip(CONDITIONS, sources, strict=True):
        with source.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != expected:
            raise ValueError(f"{source.name}: expected {expected} module pairs, found {len(rows)}")

        records = {}
        saved_pairs = {}
        points = []
        seen = set()
        for row in rows:
            record_id = row["rna_record_id"]
            if record_id not in saved_pairs:
                entry = overlap[record_id]
                names_for_pairs = entry["pair_columns"]
                rna, protein, score = [names_for_pairs.index(field) for field in (
                    "transcript_module", "proteomics_module", "weighted_jaccard",
                )]
                saved_pairs[record_id] = {(p[rna], p[protein]): p[score] for p in entry["pairs"]}
                records[record_id] = [len(records), row["cell_type"]]
            pair = row["rna_module"], row["protein_module"]
            key = record_id, *pair
            weighted_jaccard = float(row["kme_weighted_jaccard"])
            p_value = float(row["permutation_p"])
            q_value = float(row["bh_fdr_q_within_condition"])
            if (key in seen or row["combat_version"] != label or int(row["min_cells"]) != cells
                    or row["outlier_branch"] != "before_common_outlier_removal"
                    or (row["rna_specimen"], row["protein_specimen"]) != ("CSF", "CSF")
                    or int(row["hvg_count"]) != 9000 or int(row["permutations"]) != 10_000
                    or abs(saved_pairs[record_id][pair] - weighted_jaccard) > 1e-10
                    or not 0 < p_value <= q_value <= 1):
                raise ValueError(f"Unexpected scatter row for {key}")
            seen.add(key)
            points.append([records[record_id][0], *pair, weighted_jaccard, p_value, q_value,
                           int(row["shared_gene_count"])])
        datasets[correction] = {
            "label": label, "cells": cells, "stem": stem,
            "csv": (ASSET_DIR / source.name).as_posix(),
            "records": [[record_id, cell_type] for record_id, (_, cell_type) in records.items()],
            "points": points,
        }

    output = dashboard_root / ASSET_DIR
    output.mkdir(parents=True, exist_ok=True)
    for source in sources:
        shutil.copyfile(source, output / source.name)
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
