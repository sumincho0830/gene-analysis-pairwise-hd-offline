"""Validate all cross-omics module similarity deliverables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_cross_omics_kme_weighted_jaccard import (
    PROTEOMICS_REFERENCE_CONDITION,
    TRANSCRIPT_CONDITION,
    network_root,
    normalize_gene,
    output_root,
)


def proper_modules(table: pd.DataFrame) -> set[str]:
    mask = (
        ~table["is_raw_unassigned"].astype(bool)
        & table["official_assigned_abs_kME"].notna()
    )
    return set(table.loc[mask, "module_color"].astype(str))


def direct_weighted_jaccard(
    transcript: pd.DataFrame,
    protein: pd.DataFrame,
    common: set[str],
    transcript_module: str,
    protein_module: str,
) -> float:
    left = (
        transcript.loc[
            ~transcript["is_raw_unassigned"].astype(bool)
            & transcript["module_color"].astype(str).eq(transcript_module)
            & transcript["gene_symbol"].isin(common),
            ["gene_symbol", "official_assigned_abs_kME"],
        ]
        .groupby("gene_symbol")["official_assigned_abs_kME"]
        .max()
        .to_dict()
    )
    right = (
        protein.loc[
            protein["module_color"].astype(str).eq(protein_module)
            & protein["gene_symbol"].isin(common),
            ["gene_symbol", "official_assigned_abs_kME"],
        ]
        .groupby("gene_symbol")["official_assigned_abs_kME"]
        .max()
        .to_dict()
    )
    genes = set(left) | set(right)
    numerator = sum(min(float(left.get(g, 0)), float(right.get(g, 0))) for g in genes)
    denominator = sum(max(float(left.get(g, 0)), float(right.get(g, 0))) for g in genes)
    return numerator / denominator if denominator else np.nan


def validate(run_root: Path) -> dict[str, object]:
    transcript_stage = run_root / "experiments/04_variance_filtered/transcriptomics"
    cross_root = run_root / "experiments/05_cross_omics_module_similarity"
    registry = pd.read_csv(transcript_stage / "find_modules_summary.csv")
    registry = registry.loc[registry["status"].eq("complete_findModules")].copy()
    summary = pd.read_csv(cross_root / "cross_omics_similarity_summary.csv")
    errors: list[str] = []
    key_columns = ["threshold", "strategy", "mode", "specimen", "cell_type"]

    if len(summary) != len(registry):
        errors.append(f"summary records {len(summary)} != completed networks {len(registry)}")
    if summary.duplicated(key_columns).any():
        errors.append("duplicate compound experiment key in cross-omics summary")
    summary_keys = set(summary[key_columns].astype(str).itertuples(index=False, name=None))
    registry_keys = set(registry[key_columns].astype(str).itertuples(index=False, name=None))
    if summary_keys != registry_keys:
        errors.append("cross-omics compound experiment keys differ from completed registry")

    protein_tables: dict[str, pd.DataFrame] = {}
    protein_module_sets: dict[str, set[str]] = {}
    protein_universes: dict[str, set[str]] = {}
    for specimen in ("CSF", "PB"):
        source = (
            run_root
            / "experiments/01_batch_corrected/proteomics/06_kme_module_diagnosis"
            / specimen
            / PROTEOMICS_REFERENCE_CONDITION
            / "feature_module_kme.csv.gz"
        )
        original = pd.read_csv(source)
        protein_module_sets[specimen] = proper_modules(original)
        protein_universes[specimen] = set(
            original["SYMBOL"].map(normalize_gene).dropna().astype(str)
        )
        protein_tables[specimen] = pd.read_csv(
            cross_root
            / "00_shared/proteomics_gene_collapse"
            / f"{specimen}_representative_gene_modules.csv.gz"
        )

    gallery_text = (cross_root / "kme_weighted_jaccard_gallery.html").read_text(
        encoding="utf-8"
    )
    checked_direct_values = 0
    pair_rows = 0
    contribution_rows = 0
    finite_values: list[float] = []
    zero_common_modules = 0
    for row in registry.itertuples(index=False):
        matching = summary.loc[
            summary["threshold"].eq(row.threshold)
            & summary["strategy"].eq(row.strategy)
            & summary["mode"].eq(row.mode)
            & summary["specimen"].eq(row.specimen)
            & summary["cell_type"].eq(row.cell_type)
        ]
        if len(matching) != 1:
            errors.append(f"{row.network_id}: missing or duplicate summary row")
            continue
        result = matching.iloc[0]
        out = output_root(cross_root, row)
        expected_files = [
            out / "module_pair_similarity.csv",
            out / "transcript_module_best_matches.csv",
            out / "all_pair_gene_contributions.csv.gz",
            out / "reciprocal_best_gene_contributions.csv.gz",
            out / "gene_mapping_audit.json",
            out / "kme_weighted_jaccard_heatmap.png",
            out / "analysis_summary.json",
            out / "run_manifest.json",
        ]
        for path in expected_files:
            if not path.is_file() or path.stat().st_size == 0:
                errors.append(f"{row.network_id}: missing/empty {path.name}")

        pair_path = out / "module_pair_similarity.csv"
        if not pair_path.is_file():
            continue
        pairs = pd.read_csv(pair_path)
        pair_rows += len(pairs)
        all_contributions = pd.read_csv(out / "all_pair_gene_contributions.csv.gz")
        contribution_rows += len(all_contributions)
        required_contribution_columns = {
            "transcript_module",
            "proteomics_module",
            "gene",
            "transcript_abs_kME",
            "proteomics_abs_kME",
            "weighted_intersection",
            "weighted_union",
        }
        if not required_contribution_columns.issubset(all_contributions.columns):
            errors.append(f"{row.network_id}: all-pair contribution schema mismatch")
        else:
            if all_contributions.duplicated(
                ["transcript_module", "proteomics_module", "gene"]
            ).any():
                errors.append(f"{row.network_id}: duplicate pair/gene contribution")
            transcript_weight = pd.to_numeric(
                all_contributions["transcript_abs_kME"], errors="coerce"
            )
            protein_weight = pd.to_numeric(
                all_contributions["proteomics_abs_kME"], errors="coerce"
            )
            observed_min = pd.to_numeric(
                all_contributions["weighted_intersection"], errors="coerce"
            )
            observed_max = pd.to_numeric(
                all_contributions["weighted_union"], errors="coerce"
            )
            if not np.allclose(
                observed_min,
                np.minimum(transcript_weight, protein_weight),
                rtol=1e-12,
                atol=1e-12,
            ):
                errors.append(f"{row.network_id}: contribution minimum mismatch")
            if not np.allclose(
                observed_max,
                np.maximum(transcript_weight, protein_weight),
                rtol=1e-12,
                atol=1e-12,
            ):
                errors.append(f"{row.network_id}: contribution maximum mismatch")
            grouped = (
                all_contributions.groupby(
                    ["transcript_module", "proteomics_module"], sort=False
                )
                .agg(
                    contribution_shared_genes=("gene", "size"),
                    contribution_weighted_intersection=(
                        "weighted_intersection",
                        "sum",
                    ),
                )
                .reset_index()
            )
            audited = pairs.merge(
                grouped,
                on=["transcript_module", "proteomics_module"],
                how="left",
            )
            audited["contribution_shared_genes"] = audited[
                "contribution_shared_genes"
            ].fillna(0)
            audited["contribution_weighted_intersection"] = audited[
                "contribution_weighted_intersection"
            ].fillna(0.0)
            if not np.array_equal(
                audited["shared_genes"].astype(int).to_numpy(),
                audited["contribution_shared_genes"].astype(int).to_numpy(),
            ):
                errors.append(f"{row.network_id}: shared-gene counts mismatch")
            if not np.allclose(
                pd.to_numeric(audited["weighted_intersection"], errors="coerce"),
                audited["contribution_weighted_intersection"],
                rtol=1e-12,
                atol=1e-12,
            ):
                errors.append(f"{row.network_id}: pair intersection sum mismatch")
        if int(result.get("all_pair_contribution_rows", -1)) != len(all_contributions):
            errors.append(f"{row.network_id}: summary contribution row count mismatch")
        values = pd.to_numeric(pairs["kme_weighted_jaccard"], errors="coerce")
        finite = values.dropna()
        finite_values.extend(finite.astype(float).tolist())
        if ((finite < 0) | (finite > 1)).any():
            errors.append(f"{row.network_id}: similarity outside [0, 1]")

        transcript_path = (
            network_root(transcript_stage, row)
            / "06_kme_module_diagnosis"
            / TRANSCRIPT_CONDITION
            / "feature_module_kme.csv.gz"
        )
        transcript = pd.read_csv(transcript_path)
        transcript["gene_symbol"] = transcript["gene_id"].map(normalize_gene)
        expected_transcript_modules = proper_modules(transcript)
        expected_protein_modules = protein_module_sets[row.specimen]
        actual_transcript_modules = set(pairs["transcript_module"].astype(str))
        actual_protein_modules = set(pairs["proteomics_module"].astype(str))
        if actual_transcript_modules != expected_transcript_modules:
            errors.append(f"{row.network_id}: transcript module set mismatch")
        if actual_protein_modules != expected_protein_modules:
            errors.append(f"{row.network_id}: proteomics module set mismatch")
        expected_pairs = len(expected_transcript_modules) * len(expected_protein_modules)
        if len(pairs) != expected_pairs:
            errors.append(f"{row.network_id}: {len(pairs)} pairs != {expected_pairs}")

        transcript_universe = set(transcript["gene_symbol"].dropna().astype(str))
        protein = protein_tables[row.specimen]
        common = transcript_universe & protein_universes[row.specimen]
        if int(result["common_gene_universe"]) != len(common):
            errors.append(f"{row.network_id}: common universe mismatch")
        zero_common_modules += int(
            (pairs[["transcript_module", "transcript_module_genes_in_common"]]
             .drop_duplicates()["transcript_module_genes_in_common"] == 0).sum()
        )
        zero_common_modules += int(
            (pairs[["proteomics_module", "proteomics_module_genes_in_common"]]
             .drop_duplicates()["proteomics_module_genes_in_common"] == 0).sum()
        )

        best = pairs.loc[values.idxmax()]
        recalculated = direct_weighted_jaccard(
            transcript,
            protein,
            common,
            str(best["transcript_module"]),
            str(best["proteomics_module"]),
        )
        if not np.isclose(
            float(best["kme_weighted_jaccard"]), recalculated, rtol=1e-12, atol=1e-12
        ):
            errors.append(f"{row.network_id}: direct weighted Jaccard mismatch")
        checked_direct_values += 1
        figure_rel = Path(str(result["figure"])).as_posix()
        if figure_rel not in gallery_text:
            errors.append(f"{row.network_id}: figure missing from gallery")

    global_max = max(finite_values) if finite_values else np.nan
    reported_vmax = pd.to_numeric(summary["global_figure_vmax"], errors="coerce")
    if reported_vmax.nunique(dropna=False) != 1 or not np.isclose(
        float(reported_vmax.iloc[0]), global_max
    ):
        errors.append("global figure color maximum is inconsistent")

    result = {
        "status": "pass" if not errors else "fail",
        "completed_networks_expected": int(len(registry)),
        "summary_records": int(len(summary)),
        "pair_rows": int(pair_rows),
        "all_pair_gene_contribution_rows": int(contribution_rows),
        "direct_maximum_pairs_recalculated": int(checked_direct_values),
        "zero_common_module_axes_retained": int(zero_common_modules),
        "global_similarity_maximum": float(global_max),
        "proteomics_module_counts": {
            specimen: len(modules) for specimen, modules in protein_module_sets.items()
        },
        "errors": errors,
    }
    output = cross_root / "validation_summary.json"
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    args = parser.parse_args()
    validation = validate(args.run_root.resolve())
    print(json.dumps(validation, indent=2, ensure_ascii=False))
    raise SystemExit(0 if validation["status"] == "pass" else 1)
