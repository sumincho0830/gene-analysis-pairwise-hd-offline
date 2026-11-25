"""Generate transcriptomics raw QC and corrected cut=inf diagnostics.

Implemented scope:

* raw: cell-separated and cell-integrated descriptive/distribution QC;
* corrected threshold_5_cells and threshold_10_cells:
  feature provenance, CIS exclusion, before-preprocess QC, and PyWGCNA
  ``preprocess(cut=inf)`` diagnostics only.

No finite dendrogram cut, findModules, imputation, log transformation, gene
ranking, or top-N feature selection is performed.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import gc
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist, squareform

from pipeline_provenance import write_run_manifest
from transcriptomics_qc import (
    build_source_cohort_map,
    make_sample_metadata,
    read_expression_matrix,
    robust_z,
    write_expression_matrix,
    write_json,
    write_qc_outputs,
)


THRESHOLDS = ("threshold_5_cells", "threshold_10_cells")
SPECIMENS = ("CSF", "PB")


def _good_genes_fun_pandas_compat(
    datExpr: pd.DataFrame,
    weights: pd.DataFrame | None = None,
    useSamples: np.ndarray | None = None,
    useGenes: np.ndarray | None = None,
    minFraction: float = 1 / 2,
    minNSamples: int = 4,
    minNGenes: int = 4,
    tol: float | None = None,
    minRelativeWeight: float = 0.1,
) -> np.ndarray:
    """PyWGCNA 2.2.1 ``goodGenesFun`` with one pandas compatibility fix.

    PyWGCNA reconstructs ``nPresent`` in a boolean Series during the second
    filtering iteration, then assigns integer counts into it.  pandas >= 2.2
    rejects that lossy assignment.  This implementation preserves PyWGCNA's
    defaults and filtering equations, while reconstructing the counts directly
    as an integer Series.  No threshold or inclusion rule is changed.
    """
    from PyWGCNA import WGCNA

    if not datExpr.apply(lambda s: pd.to_numeric(s, errors="coerce").notnull().all()).all():
        raise ValueError("datExpr must contain numeric data.")

    weights = WGCNA.checkAndScaleWeights(weights, datExpr, scaleByMax=True)
    if tol is None:
        tol = 1e-10 * datExpr.abs().max().max()
    if useGenes is None:
        useGenes = np.repeat(True, datExpr.shape[0])
    if useSamples is None:
        useSamples = np.repeat(True, datExpr.shape[1])
    useGenes = np.asarray(useGenes, dtype=bool)
    useSamples = np.asarray(useSamples, dtype=bool)

    if len(useGenes) != datExpr.shape[0]:
        raise ValueError("Length of nGenes is not compatible with number of rows in datExpr.")
    if len(useSamples) != datExpr.shape[1]:
        raise ValueError("Length of nSamples is not compatible with number of columns in datExpr.")

    n_samples = int(useSamples.sum())
    n_genes = int(useGenes.sum())
    if weights is None:
        observed_counts = datExpr.loc[useGenes, useSamples].notna().sum(axis=1)
    else:
        observed_counts = (
            datExpr.loc[useGenes, useSamples].notna()
            & (weights.loc[useGenes, useSamples] > minRelativeWeight)
        ).sum(axis=1)

    # Compatibility fix: retain integer counts and align them to every gene.
    n_present = pd.Series(0, index=datExpr.index, dtype="int64")
    n_present.loc[observed_counts.index] = observed_counts.astype("int64")

    good_genes = useGenes.copy()
    good_genes[np.logical_and(useGenes, n_present.to_numpy() < minNSamples)] = False
    if weights is None:
        variance = np.var(datExpr.loc[good_genes, useSamples], axis=1)
    else:
        # PyWGCNA 2.2.1 does not implement weighted variance here either.
        raise NotImplementedError("Weighted goodGenesFun is not supported by PyWGCNA 2.2.1")
    variance[pd.isnull(variance)] = 0
    missing_counts = datExpr.loc[good_genes, useSamples].isna().sum(axis=1)
    good_genes[good_genes] = np.logical_and(
        np.logical_and(
            missing_counts < (1 - minFraction) * n_samples,
            variance > tol ** 2,
        ),
        n_samples - missing_counts >= minNSamples,
    )
    if int(good_genes.sum()) < minNGenes:
        raise ValueError("Too few genes with valid expression levels in the required number of samples.")
    if n_genes - int(good_genes.sum()) > 0:
        print(
            f"\n\n  ..Excluding {n_genes - int(good_genes.sum())} genes from the calculation "
            "due to too many missing samples or zero variance.\n",
            flush=True,
        )
    return good_genes


@contextmanager
def pywgcna_pandas_compatibility():
    """Apply and then restore the PyWGCNA 2.2.1 pandas compatibility shim."""
    from PyWGCNA import WGCNA

    original = WGCNA.goodGenesFun
    WGCNA.goodGenesFun = staticmethod(_good_genes_fun_pandas_compat)
    try:
        yield
    finally:
        WGCNA.goodGenesFun = staticmethod(original)


@dataclass(frozen=True)
class PipelineConfig:
    run_root: Path
    diagnostic_distance_robust_z_threshold: float = 5.0

    @property
    def formatted_root(self) -> Path:
        return self.run_root / "formatted"

    @property
    def experiment_root(self) -> Path:
        return self.run_root / "experiments"

    @property
    def clinical_path(self) -> Path:
        return self.run_root / "data/metadata/Transcriptomics_metadata.csv"


def cell_type_from_path(path: Path) -> str:
    return path.stem.removesuffix("_expression_matrix")


def integrated_paths(config: PipelineConfig, variant: str) -> tuple[Path, Path]:
    if variant == "00_raw":
        root = config.formatted_root / "00_raw/transcriptomics/pywgcna_input_cell_integrated/ALL"
        stem = "TRANSCRIPTOMICS_ALL_CELL_INTEGRATED"
    else:
        root = config.formatted_root / f"01_batch_corrected/transcriptomics/{variant}/pywgcna_input_cell_integrated/ALL"
        stem = f"TRANSCRIPTOMICS_{variant}_ALL_CELL_INTEGRATED"
    return root / f"{stem}_expression_matrix.csv", root / f"{stem}_sample_metadata.csv"


def split_integrated_formatted(config: PipelineConfig, variant: str) -> list[Path]:
    print(f"[formatted split] {variant}", flush=True)
    expression_path, metadata_path = integrated_paths(config, variant)
    matrix = read_expression_matrix(expression_path)
    metadata = pd.read_csv(metadata_path)
    key = metadata.columns[0]
    metadata[key] = metadata[key].astype(str)
    metadata = metadata.set_index(key, drop=False).loc[matrix.index]
    outputs: list[Path] = []
    base = expression_path.parents[1]
    stem = expression_path.stem.removesuffix("_expression_matrix")
    for specimen in SPECIMENS:
        ids = metadata.index[metadata["specimen"].astype(str).eq(specimen)]
        subset = matrix.loc[ids]
        subset_metadata = metadata.loc[ids]
        output_root = base / specimen
        matrix_path = output_root / f"{stem.replace('_ALL_', f'_{specimen}_')}_expression_matrix.csv"
        metadata_out = output_root / f"{stem.replace('_ALL_', f'_{specimen}_')}_sample_metadata.csv"
        write_expression_matrix(subset, matrix_path, "sample_id_specimen_cell")
        output_root.mkdir(parents=True, exist_ok=True)
        subset_metadata.to_csv(metadata_out, index=False)
        outputs.extend([matrix_path, metadata_out])
    del matrix
    gc.collect()
    return outputs


def write_feature_provenance(
    config: PipelineConfig,
    *,
    threshold: str,
    mode: str,
) -> list[Path]:
    output_root = (
        config.experiment_root
        / f"01_batch_corrected/transcriptomics/{threshold}/{mode}/00_feature_provenance"
    )
    output_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, str]] = []
    additional_rows: list[dict[str, str]] = []
    input_paths: list[Path] = []

    if mode == "cell_separate":
        raw_paths = sorted(
            (config.formatted_root / "00_raw/transcriptomics/pywgcna_input/ALL").glob("*_expression_matrix.csv")
        )
        jobs = [(cell_type_from_path(path), path) for path in raw_paths]
    else:
        raw_expression, _ = integrated_paths(config, "00_raw")
        jobs = [("cell_integrated", raw_expression)]

    for network, raw_path in jobs:
        if mode == "cell_separate":
            corrected_path = (
                config.formatted_root
                / f"01_batch_corrected/transcriptomics/{threshold}/pywgcna_input/ALL"
                / raw_path.name
            )
        else:
            corrected_path, _ = integrated_paths(config, threshold)
        raw_genes = pd.read_csv(raw_path, nrows=0).columns[1:].astype(str)
        corrected_genes = pd.read_csv(corrected_path, nrows=0).columns[1:].astype(str)
        raw_set = set(raw_genes)
        corrected_set = set(corrected_genes)
        missing = [gene for gene in raw_genes if gene not in corrected_set]
        additional = [gene for gene in corrected_genes if gene not in raw_set]
        rows.append({
            "network": network,
            "raw_gene_count": len(raw_genes),
            "corrected_gene_count": len(corrected_genes),
            "missing_from_corrected_count": len(missing),
            "additional_in_corrected_count": len(additional),
            "same_gene_order": list(raw_genes) == list(corrected_genes),
        })
        missing_rows.extend({"network": network, "gene_symbol": gene} for gene in missing)
        additional_rows.extend({"network": network, "gene_symbol": gene} for gene in additional)
        input_paths.extend([raw_path, corrected_path])

    summary_path = output_root / "feature_provenance_summary.csv"
    missing_path = output_root / "missing_genes.csv"
    additional_path = output_root / "additional_genes.csv"
    pd.DataFrame(rows).to_csv(summary_path, index=False)
    pd.DataFrame(missing_rows, columns=["network", "gene_symbol"]).to_csv(missing_path, index=False)
    pd.DataFrame(additional_rows, columns=["network", "gene_symbol"]).to_csv(additional_path, index=False)
    outputs = [summary_path, missing_path, additional_path]
    write_run_manifest(
        output_root / "run_manifest.json",
        pipeline_name=f"transcriptomics_{threshold}_{mode}_feature_provenance",
        run_root=config.run_root,
        parameters={"threshold": threshold, "mode": mode},
        generation_rules={
            "comparison": "raw and received corrected gene identifiers compared without value changes",
            "missing_gene": "present in matched raw matrix and absent from corrected matrix",
        },
        input_paths=input_paths,
        output_paths=outputs,
        code_paths=[Path(__file__)],
    )
    return [*outputs, output_root / "run_manifest.json"]


def write_qc_manifest(
    config: PipelineConfig,
    *,
    root: Path,
    pipeline_name: str,
    parameters: dict[str, Any],
    input_paths: list[Path],
    output_paths: list[Path],
) -> None:
    write_run_manifest(
        root / "run_manifest.json",
        pipeline_name=pipeline_name,
        run_root=config.run_root,
        parameters=parameters,
        generation_rules={
            "expression_values": "stored values; no additional log transform or scaling",
            "pca": "complete-case nonconstant genes; mean centering only; no imputation or top-N selection",
            "outlier": "descriptive robust z or PCA-distance robust z >= 5; flagged only",
            "heatmap": "all genes in stored order; 1st-99th percentile color range",
        },
        input_paths=input_paths,
        output_paths=output_paths,
        code_paths=[Path(__file__), Path(__file__).parent / "transcriptomics_qc.py"],
        packages=("numpy", "pandas", "matplotlib", "scikit-learn", "plotly"),
    )


def run_raw_qc(
    config: PipelineConfig,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
) -> None:
    raw_all_root = config.formatted_root / "00_raw/transcriptomics/pywgcna_input/ALL"
    experiment_root = config.experiment_root / "00_raw/transcriptomics"
    for path in sorted(raw_all_root.glob("*_expression_matrix.csv")):
        cell_type = cell_type_from_path(path)
        print(f"[raw QC] cell_separate/{cell_type}", flush=True)
        matrix = read_expression_matrix(path)
        metadata = make_sample_metadata(
            matrix.index,
            clinical=clinical,
            cohort_map=cohort_map,
            fixed_cell_type=cell_type,
        )
        cell_root = experiment_root / "cell_separate" / cell_type
        outputs = write_qc_outputs(
            matrix, metadata, cell_root,
            title=f"Raw transcriptomics: {cell_type}",
        )
        all_paths = list(outputs.paths)
        for specimen in SPECIMENS:
            ids = metadata.index[metadata["specimen"].eq(specimen)]
            if len(ids) < 3:
                continue
            specimen_root = cell_root / "02_data_distribution/by_specimen" / specimen
            specimen_outputs = write_qc_outputs(
                matrix.loc[ids], metadata.loc[ids], specimen_root,
                title=f"Raw transcriptomics: {cell_type} ({specimen})",
                include_descriptive=False,
            )
            all_paths.extend(specimen_outputs.paths)
        write_qc_manifest(
            config,
            root=cell_root,
            pipeline_name="raw_transcriptomics_cell_separate_qc",
            parameters={"cell_type": cell_type, "scope": "ALL_with_by_specimen"},
            input_paths=[path, config.clinical_path],
            output_paths=all_paths,
        )
        del matrix, metadata
        gc.collect()

    expression_path, metadata_path = integrated_paths(config, "00_raw")
    print("[raw QC] cell_integrated", flush=True)
    matrix = read_expression_matrix(expression_path)
    integrated_metadata = pd.read_csv(metadata_path)
    metadata = make_sample_metadata(
        matrix.index,
        clinical=clinical,
        cohort_map=cohort_map,
        integrated_metadata=integrated_metadata,
    )
    integrated_root = experiment_root / "cell_integrated"
    outputs = write_qc_outputs(
        matrix, metadata, integrated_root,
        title="Raw transcriptomics: cell-integrated",
    )
    all_paths = list(outputs.paths)
    for specimen in SPECIMENS:
        ids = metadata.index[metadata["specimen"].eq(specimen)]
        specimen_root = integrated_root / "02_data_distribution/by_specimen" / specimen
        specimen_outputs = write_qc_outputs(
            matrix.loc[ids], metadata.loc[ids], specimen_root,
            title=f"Raw transcriptomics: cell-integrated ({specimen})",
            include_descriptive=False,
        )
        all_paths.extend(specimen_outputs.paths)
    write_qc_manifest(
        config,
        root=integrated_root,
        pipeline_name="raw_transcriptomics_cell_integrated_qc",
        parameters={"scope": "ALL_with_by_specimen", "row_identity": "sample+specimen+cell_type"},
        input_paths=[expression_path, metadata_path, config.clinical_path],
        output_paths=all_paths,
    )
    del matrix, metadata
    gc.collect()


def grouped_distance_robust_z(metrics: pd.DataFrame, integrated: bool) -> pd.Series:
    columns = ["cohort", "cell_type"] if integrated else ["cohort"]
    result = pd.Series(np.nan, index=metrics.index, dtype=float)
    for _, indices in metrics.groupby(columns, sort=False).groups.items():
        if len(indices) >= 3:
            result.loc[indices] = robust_z(
                metrics.loc[indices, "median_pairwise_euclidean_distance"]
            )
    missing = result.isna()
    if missing.any():
        result.loc[missing] = robust_z(
            metrics.loc[missing, "median_pairwise_euclidean_distance"]
        )
    return result


def plot_diagnostic_dendrogram(
    linkage_matrix: np.ndarray,
    labels: list[str],
    output_path: Path,
    title: str,
) -> list[int]:
    width = max(18.0, min(50.0, len(labels) / 8.0))
    figure, axis = plt.subplots(figsize=(width, 9), facecolor="white")
    result = dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=5,
        color_threshold=0,
        above_threshold_color="#4C78A8",
        ax=axis,
    )
    axis.set_title(title)
    axis.set_xlabel("Samples")
    axis.set_ylabel("Euclidean distance; average linkage")
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [int(value) for value in result["leaves"]]


def write_skipped_preprocess(
    root: Path,
    *,
    matrix: pd.DataFrame,
    reason: str,
) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    removed_features = root / "removed_features.csv"
    removed_samples = root / "removed_samples.csv"
    pd.DataFrame(columns=["gene_symbol", "removal_stage", "reason"]).to_csv(removed_features, index=False)
    pd.DataFrame(columns=["sample_id", "removal_stage", "reason"]).to_csv(removed_samples, index=False)
    summary = root / "preprocess_summary.json"
    write_json(summary, {
        "status": "not_run_insufficient_samples",
        "reason": reason,
        "input_samples_after_cis_exclusion": int(matrix.shape[0]),
        "input_genes": int(matrix.shape[1]),
        "cut": "inf",
        "findModules_run": False,
    })
    return [removed_features, removed_samples, summary]


def run_cut_diagnostic(
    config: PipelineConfig,
    *,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    output_root: Path,
    network_id: str,
    input_path: Path,
    integrated: bool,
) -> list[Path]:
    from PyWGCNA import WGCNA

    diagnostic_root = output_root / "01_cut_diagnostic"
    diagnostic_root.mkdir(parents=True, exist_ok=True)
    if matrix.shape[0] < 4:
        outputs = write_skipped_preprocess(
            diagnostic_root,
            matrix=matrix,
            reason="PyWGCNA goodSamplesGenes default requires at least 4 samples",
        )
        write_run_manifest(
            diagnostic_root / "run_manifest.json",
            pipeline_name="batch_corrected_transcriptomics_cut_inf_diagnostic_skipped",
            run_root=config.run_root,
            parameters={"network_id": network_id, "cut": "inf"},
            generation_rules={"status": "explicitly skipped; no parameters relaxed"},
            input_paths=[input_path],
            output_paths=outputs,
            code_paths=[Path(__file__)],
        )
        return [*outputs, diagnostic_root / "run_manifest.json"]

    model_name = f"TRANSCRIPTOMICS_{network_id}_cut_diagnostic".replace(" ", "_")
    sample_info = metadata.loc[matrix.index, [
        "project_id", "patient_id", "cohort", "specimen", "cell_type",
        "diagnosis", "sex", "age",
    ]].copy()
    model = WGCNA(
        name=model_name,
        geneExp=matrix,
        sampleInfo=sample_info,
        species="Human",
        save=True,
        outputPath=diagnostic_root.resolve().as_posix() + "/",
        figureType="png",
        cut=float("inf"),
    )
    tpm_keep = (matrix > model.TPMcutoff).any(axis=0).to_numpy()
    after_tpm = matrix.loc[:, tpm_keep]
    with pywgcna_pandas_compatibility():
        good_genes, good_samples, _ = WGCNA.goodSamplesGenes(after_tpm.T)
        good_genes = np.asarray(good_genes, dtype=bool)
        good_samples = np.asarray(good_samples, dtype=bool)
        expected = after_tpm.loc[good_samples, good_genes].copy()
    if expected.shape[0] < 2 or expected.shape[1] == 0:
        outputs = write_skipped_preprocess(
            diagnostic_root,
            matrix=matrix,
            reason="PyWGCNA default TPMcutoff/goodSamplesGenes left an unusable matrix",
        )
        write_run_manifest(
            diagnostic_root / "run_manifest.json",
            pipeline_name="batch_corrected_transcriptomics_cut_inf_diagnostic_skipped",
            run_root=config.run_root,
            parameters={"network_id": network_id, "cut": "inf"},
            generation_rules={"status": "default preprocessing left an unusable matrix; no parameters relaxed"},
            input_paths=[input_path],
            output_paths=outputs,
            code_paths=[Path(__file__)],
        )
        return [*outputs, diagnostic_root / "run_manifest.json"]

    with pywgcna_pandas_compatibility():
        model.preprocess(show=False)
    actual = model.datExpr.to_df().copy()
    actual.index = actual.index.astype(str)
    actual.columns = actual.columns.astype(str)
    if list(actual.index) != list(expected.index) or list(actual.columns) != list(expected.columns):
        raise ValueError(f"{network_id}: PyWGCNA pre-cut IDs differ from audited expected IDs")
    if not np.array_equal(actual.to_numpy(dtype=float), expected.to_numpy(dtype=float)):
        raise ValueError(f"{network_id}: PyWGCNA pre-cut values differ from audited matrix")

    model.analysis_provenance = {
        "stage": "01_cut_diagnostic",
        "source_matrix": str(input_path),
        "cis_excluded": True,
        "transform": "none",
        "TPMcutoff": model.TPMcutoff,
        "goodSamplesGenes": (
            "PyWGCNA 2.2.1 default criteria with pandas compatibility fix; "
            "no threshold or inclusion rule changed"
        ),
        "cut": "inf",
        "findModules_run": False,
    }
    model.saveWGCNA()
    object_path = diagnostic_root / f"{model_name}.p"
    matrix_path = diagnostic_root / "matrix_before_cut.csv"
    write_expression_matrix(actual, matrix_path, "sample_id")

    removed_feature_rows = [
        {"gene_symbol": gene, "removal_stage": "TPMcutoff", "reason": f"no sample value greater than TPMcutoff={model.TPMcutoff}"}
        for gene in matrix.columns[~tpm_keep]
    ]
    removed_feature_rows.extend(
        {
            "gene_symbol": gene,
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        }
        for gene in after_tpm.columns[~good_genes]
    )
    removed_features = pd.DataFrame(
        removed_feature_rows,
        columns=["gene_symbol", "removal_stage", "reason"],
    )
    removed_features.to_csv(diagnostic_root / "removed_features.csv", index=False)
    pd.DataFrame({
        "sample_id": after_tpm.index[~good_samples].astype(str),
        "removal_stage": "goodSamplesGenes",
        "reason": "failed PyWGCNA 2.2.1 default criteria",
    }).to_csv(diagnostic_root / "removed_samples.csv", index=False)

    condensed = pdist(actual.to_numpy(dtype=float), metric="euclidean")
    linkage_matrix = WGCNA.hclust(condensed, method="average")
    pd.DataFrame(
        linkage_matrix,
        columns=["left_node", "right_node", "distance", "cluster_size"],
    ).rename_axis("merge_index").reset_index().to_csv(
        diagnostic_root / "sample_linkage.csv", index=False
    )
    leaf_order = plot_diagnostic_dendrogram(
        linkage_matrix,
        actual.index.astype(str).tolist(),
        diagnostic_root / "dendrogram_cut_diagnostic.png",
        f"Transcriptomics sample clustering ({network_id}; cut=inf diagnostic)",
    )
    pd.DataFrame({
        "leaf_position_1based": np.arange(1, len(leaf_order) + 1),
        "sample_id": actual.index[leaf_order].astype(str),
    }).to_csv(diagnostic_root / "dendrogram_leaf_order.csv", index=False)

    distances = squareform(condensed)
    np.fill_diagonal(distances, np.nan)
    metrics = metadata.loc[actual.index].copy()
    metrics["mean_pairwise_euclidean_distance"] = np.nanmean(distances, axis=1)
    metrics["median_pairwise_euclidean_distance"] = np.nanmedian(distances, axis=1)
    metrics["max_pairwise_euclidean_distance"] = np.nanmax(distances, axis=1)
    metrics["nearest_neighbor_euclidean_distance"] = np.nanmin(distances, axis=1)
    metrics["median_distance_group_robust_z"] = grouped_distance_robust_z(metrics, integrated)
    threshold = config.diagnostic_distance_robust_z_threshold
    metrics["diagnostic_distance_threshold"] = threshold
    metrics["diagnostic_outlier_candidate"] = metrics["median_distance_group_robust_z"].ge(threshold)
    metrics["candidate_action"] = np.where(metrics["diagnostic_outlier_candidate"], "flagged_only", "none")
    metrics.reset_index(drop=True).to_csv(diagnostic_root / "sample_distance_metrics.csv", index=False)
    metrics.loc[metrics["diagnostic_outlier_candidate"]].reset_index(drop=True).to_csv(
        diagnostic_root / "outlier_candidates_before_cut.csv", index=False
    )

    summary_path = diagnostic_root / "preprocess_summary.json"
    write_json(summary_path, {
        "status": "complete_cut_diagnostic",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "network_id": network_id,
        "diagnostic_only": True,
        "selected_finite_cut": None,
        "findModules_run": False,
        "input_samples_after_cis_exclusion": int(matrix.shape[0]),
        "input_genes": int(matrix.shape[1]),
        "TPMcutoff": float(model.TPMcutoff),
        "genes_removed_by_TPMcutoff": int((~tpm_keep).sum()),
        "genes_removed_by_goodSamplesGenes": int((~good_genes).sum()),
        "samples_removed_by_goodSamplesGenes": int((~good_samples).sum()),
        "cut": "inf",
        "samples_removed_by_dendrogram_cut": 0,
        "matrix_before_cut_samples": int(actual.shape[0]),
        "matrix_before_cut_genes": int(actual.shape[1]),
        "diagnostic_outlier_candidates": int(metrics["diagnostic_outlier_candidate"].sum()),
        "diagnostic_candidate_action": "flagged_only; no removal",
        "object": str(object_path),
    })
    outputs = [
        matrix_path,
        diagnostic_root / "removed_features.csv",
        diagnostic_root / "removed_samples.csv",
        diagnostic_root / "sample_linkage.csv",
        diagnostic_root / "dendrogram_leaf_order.csv",
        diagnostic_root / "dendrogram_cut_diagnostic.png",
        diagnostic_root / "sample_distance_metrics.csv",
        diagnostic_root / "outlier_candidates_before_cut.csv",
        summary_path,
        object_path,
    ]
    write_run_manifest(
        diagnostic_root / "run_manifest.json",
        pipeline_name="batch_corrected_transcriptomics_cut_inf_diagnostic",
        run_root=config.run_root,
        parameters={
            "network_id": network_id,
            "TPMcutoff": float(model.TPMcutoff),
            "cut": "inf",
            "goodSamplesGenes": (
                "PyWGCNA 2.2.1 default criteria with pandas compatibility fix; "
                "no threshold or inclusion rule changed"
            ),
        },
        generation_rules={
            "input": "CIS-excluded received corrected values",
            "transform": "none",
            "stop_gate": "diagnostic only; no finite cut and no findModules",
            "matrix_checkpoint": "TPMcutoff and goodSamplesGenes output immediately before cut",
        },
        input_paths=[input_path],
        output_paths=outputs,
        code_paths=[Path(__file__)],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )
    plt.close("all")
    del model
    gc.collect()
    return [*outputs, diagnostic_root / "run_manifest.json"]


def process_corrected_mode(
    config: PipelineConfig,
    *,
    threshold: str,
    mode: str,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
) -> None:
    mode_root = config.experiment_root / f"01_batch_corrected/transcriptomics/{threshold}/{mode}"
    write_feature_provenance(config, threshold=threshold, mode=mode)
    selection_rows: list[pd.DataFrame] = []
    selection_outputs: list[Path] = []

    if mode == "cell_separate":
        formatted_base = config.formatted_root / f"01_batch_corrected/transcriptomics/{threshold}/pywgcna_input"
        networks = []
        for specimen in SPECIMENS:
            for path in sorted((formatted_base / specimen).glob("*_expression_matrix.csv")):
                networks.append((specimen, cell_type_from_path(path), path, None))
    else:
        expression_all, metadata_all = integrated_paths(config, threshold)
        base = expression_all.parents[1]
        all_meta = pd.read_csv(metadata_all)
        networks = []
        for specimen in SPECIMENS:
            stem = expression_all.stem.removesuffix("_expression_matrix").replace("_ALL_", f"_{specimen}_")
            networks.append((
                specimen,
                "cell_integrated",
                base / specimen / f"{stem}_expression_matrix.csv",
                base / specimen / f"{stem}_sample_metadata.csv",
            ))

    for specimen, cell_type, source_path, integrated_meta_path in networks:
        print(f"[corrected] {threshold}/{mode}/{specimen}/{cell_type}", flush=True)
        matrix = read_expression_matrix(source_path)
        integrated_metadata = pd.read_csv(integrated_meta_path) if integrated_meta_path else None
        metadata = make_sample_metadata(
            matrix.index,
            clinical=clinical,
            cohort_map=cohort_map,
            fixed_cell_type=None if mode == "cell_integrated" else cell_type,
            integrated_metadata=integrated_metadata,
        )
        selection = metadata.copy()
        selection["threshold"] = threshold
        selection["mode"] = mode
        selection["network_cell_type"] = cell_type
        selection["included"] = ~selection["diagnosis"].astype(str).eq("CIS")
        selection["exclusion_reason"] = np.where(
            selection["included"], "", "diagnosis_CIS_excluded_by_design"
        )
        selection_rows.append(selection.reset_index(drop=True))
        retained_ids = selection.index[selection["included"]]
        selected = matrix.loc[retained_ids].copy()
        selected_metadata = metadata.loc[retained_ids].copy()

        matrix_name = (
            f"{cell_type}_after_cis_exclusion.csv"
            if mode == "cell_separate"
            else f"TRANSCRIPTOMICS_{threshold}_{specimen}_CELL_INTEGRATED_after_cis_exclusion.csv"
        )
        checkpoint = mode_root / "01_sample_selection/matrices" / specimen
        if mode == "cell_separate":
            checkpoint = checkpoint / cell_type
        checkpoint_path = checkpoint / matrix_name
        write_expression_matrix(
            selected,
            checkpoint_path,
            "sample_id_specimen_cell" if mode == "cell_integrated" else "gene_symbol",
        )
        selection_outputs.append(checkpoint_path)

        before_root = mode_root / "02_before_preprocess" / specimen
        preprocess_root = mode_root / "03_preprocess" / specimen
        if mode == "cell_separate":
            before_root = before_root / cell_type
            preprocess_root = preprocess_root / cell_type
        qc_manifest_path = before_root / "run_manifest.json"
        preprocess_manifest_path = preprocess_root / "01_cut_diagnostic/run_manifest.json"
        if qc_manifest_path.exists() and preprocess_manifest_path.exists():
            print(
                f"[resume skip] completed {threshold}/{mode}/{specimen}/{cell_type}",
                flush=True,
            )
            del matrix, selected, metadata, selected_metadata
            gc.collect()
            continue
        qc_outputs = write_qc_outputs(
            selected,
            selected_metadata,
            before_root,
            title=f"Corrected transcriptomics {threshold}: {cell_type} ({specimen})",
        )
        write_qc_manifest(
            config,
            root=before_root,
            pipeline_name=f"batch_corrected_transcriptomics_{mode}_before_preprocess_qc",
            parameters={
                "threshold": threshold,
                "mode": mode,
                "specimen": specimen,
                "cell_type": cell_type,
                "cis_excluded": True,
            },
            input_paths=[checkpoint_path, config.clinical_path],
            output_paths=list(qc_outputs.paths),
        )
        run_cut_diagnostic(
            config,
            matrix=selected,
            metadata=selected_metadata,
            output_root=preprocess_root,
            network_id=f"{threshold}_{mode}_{specimen}_{cell_type}",
            input_path=checkpoint_path,
            integrated=mode == "cell_integrated",
        )
        del matrix, selected, metadata, selected_metadata
        gc.collect()

    selection_root = mode_root / "01_sample_selection"
    selection_root.mkdir(parents=True, exist_ok=True)
    selection_table = pd.concat(selection_rows, ignore_index=True)
    selection_path = selection_root / "sample_selection.csv"
    excluded_path = selection_root / "excluded_cis_samples.csv"
    retained_path = selection_root / "retained_samples.csv"
    counts_path = selection_root / "network_sample_counts.csv"
    selection_table.to_csv(selection_path, index=False)
    selection_table.loc[~selection_table["included"]].to_csv(excluded_path, index=False)
    selection_table.loc[selection_table["included"]].to_csv(retained_path, index=False)
    selection_table.groupby(
        ["threshold", "mode", "specimen", "network_cell_type"], dropna=False
    )["included"].agg(samples_before="size", samples_after="sum").reset_index().assign(
        cis_excluded=lambda x: x["samples_before"] - x["samples_after"]
    ).to_csv(counts_path, index=False)
    selection_outputs.extend([selection_path, excluded_path, retained_path, counts_path])
    write_run_manifest(
        selection_root / "run_manifest.json",
        pipeline_name=f"batch_corrected_transcriptomics_{threshold}_{mode}_cis_selection",
        run_root=config.run_root,
        parameters={"excluded_diagnosis": "CIS", "threshold": threshold, "mode": mode},
        generation_rules={
            "selection": "exclude rows whose source-sample clinical diagnosis is exactly CIS",
            "batch_correction_refit": "not performed; received corrected values retained",
            "checkpoint": "one CIS-excluded matrix per downstream network",
        },
        input_paths=[config.clinical_path],
        output_paths=selection_outputs,
        code_paths=[Path(__file__)],
    )


def run(
    config: PipelineConfig,
    *,
    only_threshold: str | None = None,
    skip_raw: bool = False,
    skip_formatted_splits: bool = False,
) -> None:
    clinical = pd.read_csv(config.clinical_path)
    cohort_map = build_source_cohort_map(config.run_root)
    if not skip_formatted_splits:
        generated_formatted: list[Path] = []
        variants = THRESHOLDS if skip_raw else ("00_raw", *THRESHOLDS)
        for variant in variants:
            generated_formatted.extend(split_integrated_formatted(config, variant))
        write_json(
            config.formatted_root / "manifests/transcriptomics_generated_splits.json",
            {
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "rule": "split cell-integrated ALL matrix by metadata specimen without changing values",
                "outputs": [str(path) for path in generated_formatted],
            },
        )
    if not skip_raw:
        run_raw_qc(config, clinical, cohort_map)
    selected_thresholds = (only_threshold,) if only_threshold else THRESHOLDS
    for threshold in selected_thresholds:
        for mode in ("cell_separate", "cell_integrated"):
            process_corrected_mode(
                config,
                threshold=threshold,
                mode=mode,
                clinical=clinical,
                cohort_map=cohort_map,
            )


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Generate transcriptomics raw QC and corrected cut=inf diagnostics."
    )
    parser.add_argument("--run-root", type=Path, default=default_root)
    parser.add_argument("--only-threshold", choices=THRESHOLDS)
    parser.add_argument("--skip-raw", action="store_true")
    parser.add_argument("--skip-formatted-splits", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = PipelineConfig(run_root=args.run_root.resolve())
    run(
        config,
        only_threshold=args.only_threshold,
        skip_raw=args.skip_raw,
        skip_formatted_splits=args.skip_formatted_splits,
    )
    print("Completed transcriptomics raw QC and corrected cut=inf diagnostics.")


if __name__ == "__main__":
    main()

