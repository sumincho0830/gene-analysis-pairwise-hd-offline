"""Build the offline module-pair viewer data from the shipped CSV files."""

import argparse
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]


def read_rows(record, file_id, dashboard_root):
    path = next((f["path"] for f in record["files"] if f["id"] == file_id), None)
    if path:
        with (dashboard_root / path).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    return []


def number(value):
    return float(value) if value else None


def build(dashboard_root: Path = ROOT):
    source = (dashboard_root / "dashboard-data.js").read_text(encoding="utf-8")
    records = json.loads(source.split("=", 1)[1].strip().removesuffix(";"))["records"]
    pairs = {}
    diagnosis = {}
    for record in records:
        record_id = record["id"]
        rows = read_rows(record, "cross_pair_similarity", dashboard_root)
        if rows:
            matrix = {}
            for row in rows:
                matrix.setdefault(row["transcript_module"], {})[row["proteomics_module"]] = [
                    float(row["kme_weighted_jaccard"]),
                    float(row["jaccard"]),
                    2 if row["reciprocal_best"] == "True" else int(row["transcript_best"] == "True"),
                ]
            pairs[record_id] = matrix
        diagnosis_sources = [(record_id, "diagnosis_pairwise_statistics")]
        if record["dimensions"]["mode"] != "proteomics":
            diagnosis_sources.append(
                (f"{record_id}::after_outlier", "after_outlier_diagnosis_pairwise_statistics")
            )
        for diagnosis_id, file_id in diagnosis_sources:
            for row in read_rows(record, file_id, dashboard_root):
                modules = diagnosis.setdefault(diagnosis_id, {})
                modules.setdefault(row["module_color"], []).append([
                    row["contrast"],
                    number(row["point_biserial_r"]),
                    number(row["p_value_two_sided"]),
                    number(row["bh_fdr_q_within_network"]),
                    int(row["n_reference"]),
                    int(row["n_target"]),
                    int(row["n_pairwise"]),
                    row["inference_ready"].lower() == "true",
                    row["non_estimable_reason"],
                ])
    maximum = max(
        (cell[0] for matrix in pairs.values() for row in matrix.values() for cell in row.values()),
        default=0.0,
    )
    payload = json.dumps(
        {"pairs": pairs, "diagnosis": diagnosis, "maxWeightedJaccard": maximum},
        separators=(",", ":"),
        allow_nan=False,
    )
    (dashboard_root / "module-members/cross-omics-view.js").write_bytes(
        f"window.CROSS_OMICS_VIEW = {payload};\n".encode("utf-8")
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dashboard-root", type=Path, default=ROOT)
    args = parser.parse_args()
    build(args.dashboard_root.resolve())
