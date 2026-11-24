"""Run official PyWGCNA kME and module-diagnosis analysis after findModules.

This stage reconstructs a light analysis object from the audited selected-cut
expression object plus the saved final module assignments/eigengenes.  It then
calls the installed PyWGCNA 2.2.1 public methods directly:

* ``WGCNA.CalculateSignedKME()``
* ``WGCNA.module_trait_relationships_heatmap(metaData=['diagnosis'])``

``CalculateSignedKME`` internally allocates correlations among every supplied
feature.  Calling it on independent feature chunks is mathematically identical
for the requested feature-to-eigengene block and avoids an unused O(p^2)
feature-feature allocation for 5,000--12,000-gene networks.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import gc
import gzip
import json
from pathlib import Path
import pickle
import shutil
import subprocess
from typing import Any, Iterator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pipeline_provenance import write_run_manifest


CONDITION = "01_default_ds2_mms50_med0p20"
PROTEOMICS_CONDITIONS = (
    "01_default_ds2_mms50_med0p20",
    "02_ds3_mms50_med0p20",
    "03_ds3_mms30_med0p20",
    "04_ds3_mms30_med0p10",
)
KME_CHUNK_SIZE = 1000


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


@contextmanager
def open_pickle_stream(path: Path) -> Iterator[Any]:
    if path.name.endswith(".p.gz"):
        with gzip.open(path, "rb") as stream:
            yield stream
        return
    if path.name.endswith(".p.zst"):
        zstd = shutil.which("zstd")
        if not zstd:
            raise RuntimeError("zstd executable is required to load .p.zst objects")
        with path.open("rb") as source:
            process = subprocess.Popen(
                [zstd, "-dc", "-q"], stdin=source, stdout=subprocess.PIPE
            )
            if process.stdout is None:
                raise RuntimeError(f"Could not open zstd stream: {path}")
            try:
                yield process.stdout
            finally:
                process.stdout.close()
                return_code = process.wait()
                if return_code != 0:
                    raise RuntimeError(f"zstd decompression failed: {path}")
        return
    if path.name.endswith(".p"):
        with path.open("rb") as stream:
            yield stream
        return
    raise ValueError(f"Unsupported PyWGCNA object: {path}")


def load_model(path: Path) -> Any:
    with open_pickle_stream(path) as stream:
        return pickle.load(stream)


def load_analysis_model(
    *,
    source_model: Path | None,
    source_matrix: Path | None = None,
    sample_info_path: Path | None = None,
) -> Any:
    """Load a stored model or reconstruct a light object from audited tables."""

    if source_model is not None:
        return load_model(source_model)
    if source_matrix is None or sample_info_path is None:
        raise ValueError("source_matrix and sample_info_path are required without a model")
    from PyWGCNA import WGCNA

    expression = pd.read_csv(source_matrix, index_col="sample_id")
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)
    sample_info = pd.read_csv(sample_info_path, index_col="sample_id")
    sample_info.index = sample_info.index.astype(str)
    sample_info = sample_info.reindex(expression.index)
    if sample_info.isna().all(axis=1).any():
        raise ValueError("Sample metadata does not cover the reconstructed expression matrix")
    return WGCNA(
        name="TRANSCRIPTOMICS_AFTER_OUTLIER_REMOVAL_DOWNSTREAM",
        geneExp=expression,
        sampleInfo=sample_info,
        species="Human",
        save=False,
    )


def find_model(directory: Path) -> Path:
    candidates = sorted(
        path
        for pattern in ("*.p", "*.p.gz", "*.p.zst")
        for path in directory.glob(pattern)
    )
    if len(candidates) != 1:
        raise ValueError(f"Expected one selected-cut object in {directory}: {candidates}")
    return candidates[0]


def find_model_used_by_findmodules(
    run_root: Path, module_root: Path, selected_cut_root: Path
) -> Path:
    manifest_path = module_root / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for record in manifest.get("inputs", []):
            recorded = str(record.get("path", ""))
            if recorded.endswith((".p", ".p.gz", ".p.zst")):
                candidate = run_root / Path(recorded)
                if candidate.exists():
                    return candidate
    return find_model(selected_cut_root)


def read_eigengenes(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    if "sample_id" not in table.columns:
        raise ValueError(f"sample_id is missing from {path}")
    eigengene_columns = [column for column in table.columns if column.startswith("ME")]
    if not eigengene_columns:
        raise ValueError(f"No ME columns in {path}")
    result = table.set_index("sample_id")[eigengene_columns].astype(float)
    result.index = result.index.astype(str)
    return result


def read_assignments(path: Path, feature_column: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    if feature_column not in table.columns or "module_color" not in table.columns:
        raise ValueError(f"Required assignment columns are missing from {path}")
    table[feature_column] = table[feature_column].astype(str)
    return table.set_index(feature_column, drop=False)


def prepare_model(
    model: Any,
    assignments: pd.DataFrame,
    eigengenes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    feature_ids = pd.Index(model.datExpr.var_names.astype(str), name="feature_id")
    sample_ids = pd.Index(model.datExpr.obs_names.astype(str), name="sample_id")
    if feature_ids.has_duplicates or sample_ids.has_duplicates:
        raise ValueError("Selected-cut PyWGCNA object has duplicated identifiers")

    assignments = assignments.reindex(feature_ids)
    if assignments["module_color"].isna().any():
        missing = assignments.index[assignments["module_color"].isna()].tolist()[:10]
        raise ValueError(f"Module assignments do not cover expression features: {missing}")
    eigengenes = eigengenes.reindex(sample_ids)
    if eigengenes.isna().to_numpy().any():
        raise ValueError("Module eigengenes do not cover expression samples")

    expected = {"ME" + color for color in assignments["module_color"].astype(str).unique()}
    if expected != set(eigengenes.columns):
        raise ValueError(
            f"Assignment/eigengene modules differ: expected={sorted(expected)}, "
            f"observed={sorted(eigengenes.columns)}"
        )

    raw_unassigned_colors: list[str] = []
    if "is_raw_unassigned" in assignments.columns:
        for color, group in assignments.groupby("module_color", sort=False):
            if group["is_raw_unassigned"].astype(bool).all():
                raw_unassigned_colors.append(str(color))
    excluded_columns = {"ME" + color for color in raw_unassigned_colors}
    proper_eigengenes = eigengenes.drop(
        columns=[column for column in eigengenes.columns if column in excluded_columns]
    )
    if proper_eigengenes.empty:
        raise ValueError("No proper module eigengenes remain after raw-unassigned exclusion")

    model.datExpr.var["moduleColors"] = assignments["module_color"].astype(str).to_numpy()
    model.datME = proper_eigengenes.copy()
    model.MEs = proper_eigengenes.copy()
    return assignments, eigengenes, proper_eigengenes, raw_unassigned_colors


def official_signed_kme_chunked(
    model: Any,
    eigengenes: pd.DataFrame,
    chunk_size: int = KME_CHUNK_SIZE,
) -> pd.DataFrame:
    """Call PyWGCNA.CalculateSignedKME on independent feature chunks."""

    original_dat_expr = model.datExpr
    parts: list[pd.DataFrame] = []
    try:
        model.datME = eigengenes
        model.MEs = eigengenes
        for start in range(0, original_dat_expr.n_vars, chunk_size):
            stop = min(start + chunk_size, original_dat_expr.n_vars)
            model.datExpr = original_dat_expr[:, start:stop].copy()
            model.CalculateSignedKME()
            parts.append(model.signedKME.copy())
    finally:
        model.datExpr = original_dat_expr
        model.datME = eigengenes
        model.MEs = eigengenes

    result = pd.concat(parts, axis=0)
    result.index = result.index.astype(str)
    result.index.name = "feature_id"
    result = result.reindex(original_dat_expr.var_names.astype(str))
    if result.isna().all(axis=1).any():
        raise ValueError("Official kME calculation produced all-NA feature rows")
    model.signedKME = result
    return result


def bh_fdr(values: np.ndarray) -> np.ndarray:
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    pvalues = values[finite]
    if not len(pvalues):
        return result
    order = np.argsort(pvalues, kind="stable")
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[finite] = restored
    return result


def official_module_diagnosis(
    model: Any,
    output_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """Run the public PyWGCNA module-trait heatmap method for diagnosis."""

    if "diagnosis" not in model.datExpr.obs.columns:
        raise ValueError("diagnosis is missing from PyWGCNA sample metadata")
    diagnosis = model.datExpr.obs["diagnosis"].astype(str)
    if diagnosis.nunique() < 2:
        empty = pd.DataFrame()
        encoding = pd.DataFrame(
            {"sample_id": model.datExpr.obs_names.astype(str), "diagnosis": diagnosis.to_numpy()}
        )
        return empty, empty, empty, encoding, "not_run_single_diagnosis"

    original_obs = model.datExpr.obs.copy()
    model.datExpr.obs = original_obs[["diagnosis"]].copy()
    # PyWGCNA 2.2.1 constructs the correct one-hot values but, with pandas 3,
    # leaves their dtype as object. SciPy 1.17 pearsonr rejects object arrays.
    # Casting that official encoding to float preserves every value and lets the
    # public heatmap method follow its numeric-trait branch.
    encoding_values = model.getDatTraits(["diagnosis"]).astype(float)
    encoding = pd.concat(
        [model.datExpr.obs[["diagnosis"]].copy(), encoding_values.copy()], axis=1
    ).reset_index().rename(columns={"index": "sample_id"})
    if encoding.columns[0] != "sample_id":
        encoding = encoding.rename(columns={encoding.columns[0]: "sample_id"})

    figures_root = output_root / "figures"
    figures_root.mkdir(parents=True, exist_ok=True)
    model.name = model.name or "WGCNA"
    model.outputPath = output_root.resolve().as_posix() + "/"
    model.figureType = "png"
    model.save = True
    try:
        model.datExpr.obs = encoding_values.copy()
        model.module_trait_relationships_heatmap(
            metaData=encoding_values.columns.astype(str).tolist(),
            alternative="two-sided",
            figsize=None,
            show=False,
            file_name="module_diagnosis_heatmap",
        )
        correlations = model.moduleTraitCor.copy()
        pvalues = model.moduleTraitPvalue.copy()
    finally:
        model.datExpr.obs = original_obs
        plt.close("all")

    correlations.index.name = "module_eigengene"
    pvalues.index.name = "module_eigengene"
    rows: list[dict[str, Any]] = []
    flattened_p = pvalues.to_numpy(dtype=float).reshape(-1)
    qvalues = bh_fdr(flattened_p).reshape(pvalues.shape)
    for i, module in enumerate(correlations.index):
        for j, trait in enumerate(correlations.columns):
            rows.append(
                {
                    "module_eigengene": module,
                    "module_color": str(module)[2:] if str(module).startswith("ME") else str(module),
                    "diagnosis_trait": trait,
                    "pearson_r": correlations.iat[i, j],
                    "p_value_two_sided": pvalues.iat[i, j],
                    "bh_fdr_q_within_network": qvalues[i, j],
                    "n_samples": int(model.datExpr.n_obs),
                }
            )
    return correlations, pvalues, pd.DataFrame(rows), encoding, "complete"


def build_feature_kme(
    assignments: pd.DataFrame,
    signed_kme: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata_columns = [
        column
        for column in assignments.columns
        if not column.startswith("kME")
        and "kme" not in column.lower()
        and column != "hub_rank_within_module"
    ]
    result = assignments[metadata_columns].copy()
    result.insert(0, "feature_id", result.index.astype(str))
    result = result.join(signed_kme, how="left")
    colors = assignments["module_color"].astype(str)
    official_module_kme = np.full(len(result), np.nan, dtype=float)
    for color in colors.unique():
        column = "kME" + color
        if column not in signed_kme.columns:
            color_rows = assignments.loc[colors.eq(color)]
            if (
                "is_raw_unassigned" in color_rows.columns
                and color_rows["is_raw_unassigned"].astype(bool).all()
            ):
                continue
            raise ValueError(f"Missing official kME column for proper module {color}")
        mask = colors.eq(color).to_numpy()
        official_module_kme[mask] = signed_kme.loc[mask, column].to_numpy(dtype=float)
    result["official_module_kME"] = official_module_kme

    raw_unassigned = (
        assignments["is_raw_unassigned"].astype(bool).to_numpy()
        if "is_raw_unassigned" in assignments.columns
        else np.zeros(len(assignments), dtype=bool)
    )
    result["official_assigned_signed_kME"] = np.where(
        raw_unassigned, np.nan, official_module_kme
    )
    result["official_assigned_abs_kME"] = np.abs(
        result["official_assigned_signed_kME"]
    )
    result["hub_rank_within_module"] = (
        result.groupby("module_color", sort=False)["official_assigned_abs_kME"]
        .rank(method="first", ascending=False)
        .astype("Int64")
    )

    summary_rows: list[dict[str, Any]] = []
    for color, group in result.groupby("module_color", sort=False):
        summary_rows.append(
            {
                "module_color": color,
                "feature_count": int(len(group)),
                "raw_unassigned_feature_count": int(
                    group.get("is_raw_unassigned", pd.Series(False, index=group.index))
                    .astype(bool)
                    .sum()
                ),
                "features_with_official_assigned_kME": int(
                    group["official_assigned_signed_kME"].notna().sum()
                ),
                "median_official_assigned_abs_kME": group[
                    "official_assigned_abs_kME"
                ].median(),
                "max_official_assigned_abs_kME": group[
                    "official_assigned_abs_kME"
                ].max(),
            }
        )
    return result, pd.DataFrame(summary_rows)


def output_complete(output_root: Path) -> bool:
    required = [
        output_root / "feature_module_kme.csv.gz",
        output_root / "module_kme_summary.csv",
        output_root / "module_diagnosis_correlations.csv",
        output_root / "module_diagnosis_pvalues.csv",
        output_root / "module_diagnosis_statistics.csv",
        output_root / "diagnosis_trait_encoding.csv",
        output_root / "figures" / "module_diagnosis_heatmap.png",
        output_root / "analysis_summary.json",
        output_root / "run_manifest.json",
    ]
    return all(path.exists() and path.stat().st_size > 0 for path in required)


def run_network(
    *,
    run_root: Path,
    omics: str,
    network_id: str,
    source_model: Path | None,
    assignments_path: Path,
    assignment_feature_column: str,
    eigengenes_path: Path,
    output_root: Path,
    dimensions: dict[str, str],
    overwrite: bool,
    source_matrix: Path | None = None,
    sample_info_path: Path | None = None,
) -> dict[str, Any]:
    if output_complete(output_root) and not overwrite:
        summary = json.loads((output_root / "analysis_summary.json").read_text(encoding="utf-8"))
        summary["resume_action"] = "reused_complete_output"
        return summary

    output_root.mkdir(parents=True, exist_ok=True)
    model = load_analysis_model(
        source_model=source_model,
        source_matrix=source_matrix,
        sample_info_path=sample_info_path,
    )
    assignments = read_assignments(assignments_path, assignment_feature_column)
    eigengenes = read_eigengenes(eigengenes_path)
    assignments, stored_eigengenes, eigengenes, raw_unassigned_colors = prepare_model(
        model, assignments, eigengenes
    )

    signed_kme = official_signed_kme_chunked(model, eigengenes)
    feature_kme, module_kme_summary = build_feature_kme(assignments, signed_kme)
    feature_kme_path = output_root / "feature_module_kme.csv.gz"
    feature_kme.to_csv(
        feature_kme_path,
        index=False,
        compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
    )
    module_kme_summary_path = output_root / "module_kme_summary.csv"
    module_kme_summary.to_csv(module_kme_summary_path, index=False)

    display_cell = dimensions.get("cell_type", "ALL")
    display_strategy = dimensions.get("strategy", dimensions.get("condition", "modules"))
    model.name = f"{omics.upper()}_{dimensions.get('specimen', 'ALL')}_{display_cell}_{display_strategy}"
    correlations, pvalues, statistics, encoding, heatmap_status = official_module_diagnosis(
        model, output_root
    )
    correlation_path = output_root / "module_diagnosis_correlations.csv"
    pvalue_path = output_root / "module_diagnosis_pvalues.csv"
    statistics_path = output_root / "module_diagnosis_statistics.csv"
    encoding_path = output_root / "diagnosis_trait_encoding.csv"
    correlations.to_csv(correlation_path)
    pvalues.to_csv(pvalue_path)
    statistics.to_csv(statistics_path, index=False)
    encoding.to_csv(encoding_path, index=False)

    summary = {
        "status": "complete" if heatmap_status == "complete" else "complete_kme_heatmap_not_run",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "omics": omics,
        "network_id": network_id,
        **dimensions,
        "samples": int(model.datExpr.n_obs),
        "features": int(model.datExpr.n_vars),
        "module_eigengenes": int(eigengenes.shape[1]),
        "stored_module_eigengenes": int(stored_eigengenes.shape[1]),
        "excluded_raw_unassigned_eigengenes": raw_unassigned_colors,
        "kme_values": int(signed_kme.size),
        "kme_method": "PyWGCNA.WGCNA.CalculateSignedKME",
        "kme_chunk_size": KME_CHUNK_SIZE,
        "diagnosis_levels": sorted(model.datExpr.obs["diagnosis"].astype(str).unique().tolist()),
        "diagnosis_traits": correlations.columns.astype(str).tolist(),
        "module_diagnosis_tests": int(statistics.shape[0]),
        "non_estimable_module_diagnosis_tests": int(
            statistics["p_value_two_sided"].isna().sum()
        ) if not statistics.empty else 0,
        "non_estimable_module_eigengenes": (
            statistics.loc[statistics["p_value_two_sided"].isna(), "module_eigengene"]
            .drop_duplicates()
            .astype(str)
            .tolist()
            if not statistics.empty
            else []
        ),
        "minimum_p_value": (
            float(statistics["p_value_two_sided"].min()) if not statistics.empty else None
        ),
        "minimum_bh_fdr_q_within_network": (
            float(statistics["bh_fdr_q_within_network"].min())
            if not statistics.empty
            else None
        ),
        "heatmap_method": "PyWGCNA.WGCNA.module_trait_relationships_heatmap",
        "heatmap_status": heatmap_status,
        "official_heatmap_annotations": "Pearson r and two-sided unadjusted p-value",
        "supplemental_multiple_testing": "BH-FDR across module x diagnosis tests within this network",
        "non_estimable_policy": "constant or numerically near-constant module eigengenes remain NA; never replace with zero",
    }
    summary_path = output_root / "analysis_summary.json"
    write_json(summary_path, summary)

    output_paths = [
        feature_kme_path,
        module_kme_summary_path,
        correlation_path,
        pvalue_path,
        statistics_path,
        encoding_path,
        summary_path,
    ]
    heatmap_path = output_root / "figures" / "module_diagnosis_heatmap.png"
    if heatmap_path.exists():
        output_paths.append(heatmap_path)
    write_run_manifest(
        output_root / "run_manifest.json",
        pipeline_name=f"{omics}_official_kme_module_diagnosis",
        run_root=run_root,
        parameters={
            "trait": "diagnosis",
            "alternative": "two-sided",
            "kme_chunk_size": KME_CHUNK_SIZE,
            **dimensions,
        },
        generation_rules={
            "expression": (
                "stored selected-cut PyWGCNA datExpr"
                if source_model is not None
                else "audited cut=inf matrix after Ng_MOG1 and Ng_HD3 removal"
            ) + "; no additional transform/filter/imputation",
            "module_assignment": "saved final findModules module_color aligned by feature ID",
            "module_eigengene": "saved final findModules ME matrix aligned by sample ID; raw-label-0 fully unassigned eigengene excluded like intended MEgrey handling",
            "kME": "public PyWGCNA.CalculateSignedKME called per independent feature chunk; all feature-to-ME correlations retained",
            "diagnosis_encoding": "public PyWGCNA.getDatTraits values; object dtype cast to float for pandas 3 / SciPy 1.17 compatibility",
            "heatmap": "public PyWGCNA.module_trait_relationships_heatmap with diagnosis and two-sided Pearson test",
            "heatmap_annotation": "official unadjusted Pearson p-value; supplemental BH-FDR is CSV only",
        },
        input_paths=[
            path
            for path in (
                source_model,
                source_matrix,
                sample_info_path,
                assignments_path,
                eigengenes_path,
            )
            if path is not None
        ],
        output_paths=output_paths,
        code_paths=[Path(__file__), Path(__file__).parent / "pipeline_provenance.py"],
        packages=("numpy", "pandas", "matplotlib", "scipy", "seaborn", "PyWGCNA", "anndata"),
    )
    del model, signed_kme, feature_kme, correlations, pvalues, statistics
    gc.collect()
    return summary


def run_proteomics(run_root: Path, overwrite: bool, limit: int | None) -> list[dict[str, Any]]:
    experiment = run_root / "experiments" / "01_batch_corrected" / "proteomics"
    summaries: list[dict[str, Any]] = []
    for specimen in ("CSF", "PB"):
        source_model = (
            experiment / "03_preprocess" / specimen / "02_selected_cut" / f"PROT_{specimen}_selected_cut.p"
        )
        for condition in PROTEOMICS_CONDITIONS:
            if limit is not None and len(summaries) >= limit:
                return summaries
            module_root = experiment / "05_find_modules" / specimen / condition
            output_root = experiment / "06_kme_module_diagnosis" / specimen / condition
            print(f"[proteomics] {specimen}/{condition}", flush=True)
            summaries.append(
                run_network(
                    run_root=run_root,
                    omics="proteomics",
                    network_id=f"proteomics_{specimen}_{condition}",
                    source_model=source_model,
                    assignments_path=module_root / "feature_module_kme.csv",
                    assignment_feature_column="PROBEID",
                    eigengenes_path=module_root / "module_eigengenes.csv",
                    output_root=output_root,
                    dimensions={"specimen": specimen, "condition": condition},
                    overwrite=overwrite,
                )
            )
    return summaries


def run_transcriptomics(
    run_root: Path, overwrite: bool, limit: int | None
) -> list[dict[str, Any]]:
    experiment = run_root / "experiments" / "04_variance_filtered" / "transcriptomics"
    registry = pd.read_csv(experiment / "find_modules_summary.csv")
    registry = registry.loc[registry["status"].eq("complete_findModules")].copy()
    summaries: list[dict[str, Any]] = []
    for row in registry.itertuples(index=False):
        if limit is not None and len(summaries) >= limit:
            break
        network_root = experiment / row.threshold / row.strategy / row.mode / row.specimen
        if row.mode == "cell_separate":
            network_root = network_root / row.cell_type
        module_root = network_root / "05_find_modules" / CONDITION
        selected_cut_root = network_root / "03_preprocess" / "02_selected_cut"
        output_root = network_root / "06_kme_module_diagnosis" / CONDITION
        print(f"[transcriptomics] {row.network_id} / {row.strategy}", flush=True)
        summaries.append(
            run_network(
                run_root=run_root,
                omics="transcriptomics",
                network_id=str(row.network_id),
                source_model=find_model_used_by_findmodules(
                    run_root, module_root, selected_cut_root
                ),
                assignments_path=module_root / "module_assignments.csv.gz",
                assignment_feature_column="gene_id",
                eigengenes_path=module_root / "module_eigengenes.csv",
                output_root=output_root,
                dimensions={
                    "threshold": str(row.threshold),
                    "strategy": str(row.strategy),
                    "mode": str(row.mode),
                    "specimen": str(row.specimen),
                    "cell_type": str(row.cell_type),
                    "condition": CONDITION,
                },
                overwrite=overwrite,
            )
        )
    return summaries


def write_global_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument(
        "--omics", choices=("all", "proteomics", "transcriptomics"), default="all"
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_root = args.run_root.resolve()
    if args.omics in ("all", "proteomics"):
        rows = run_proteomics(run_root, args.overwrite, args.limit)
        write_global_summary(
            run_root
            / "experiments"
            / "01_batch_corrected"
            / "proteomics"
            / "06_kme_module_diagnosis"
            / "summary.csv",
            rows,
        )
    if args.omics in ("all", "transcriptomics"):
        rows = run_transcriptomics(run_root, args.overwrite, args.limit)
        write_global_summary(
            run_root
            / "experiments"
            / "04_variance_filtered"
            / "transcriptomics"
            / "module_diagnosis_summary.csv",
            rows,
        )


if __name__ == "__main__":
    main()
