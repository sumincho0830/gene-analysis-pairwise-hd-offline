"""Apply reviewed transcriptomics sample-tree cuts and stop before findModules.

Each network is rebuilt from its CIS-excluded checkpoint with a new PyWGCNA
object.  The TPMcutoff/goodSamplesGenes output must exactly match the preceding
``01_cut_diagnostic/matrix_before_cut.csv`` checkpoint.  Reviewed ``use`` rows
apply their finite cut and retain the unique largest cluster; reviewed
``reject`` rows and networks with zero diagnostic flags use ``cut=inf``.

PyWGCNA 2.2.1 always retains cluster label 0, although scipy's label numbers
depend on input order.  For a finite cut only, this runner swaps the largest
cluster label with 0 inside ``preprocess()`` when necessary.  The selected cut,
original labels, kept label, and whether relabelling occurred are all saved.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist

from pipeline_provenance import write_run_manifest
from run_transcriptomics_stage1 import (
    PipelineConfig,
    THRESHOLDS,
    build_source_cohort_map,
    integrated_paths,
    make_sample_metadata,
    pywgcna_pandas_compatibility,
    read_expression_matrix,
    write_expression_matrix,
    write_json,
)


MODES = ("cell_separate", "cell_integrated")
SPECIMENS = ("CSF", "PB")
REVIEW_RELATIVE_PATH = Path(
    "experiments/01_batch_corrected/transcriptomics/dendrogram/selected_cuts_review.csv"
)


@dataclass(frozen=True)
class Network:
    threshold: str
    mode: str
    specimen: str
    cell_type: str
    checkpoint_path: Path
    diagnostic_root: Path
    selected_root: Path

    @property
    def key(self) -> tuple[str, str, str, str]:
        return self.threshold, self.mode, self.specimen, self.cell_type

    @property
    def network_id(self) -> str:
        return "_".join(self.key)


@dataclass(frozen=True)
class CutDecision:
    decision: str
    selected_cut: float | None
    source: str
    confirmed: bool
    notes: str

    @property
    def cut_applied(self) -> bool:
        return self.decision == "use"

    @property
    def model_cut(self) -> float:
        return float(self.selected_cut) if self.cut_applied else float("inf")

    @property
    def serialized_cut(self) -> float | str:
        return float(self.selected_cut) if self.cut_applied else "inf"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_networks(config: PipelineConfig) -> Iterator[Network]:
    transcript_root = config.experiment_root / "01_batch_corrected/transcriptomics"
    for threshold in THRESHOLDS:
        for mode in MODES:
            mode_root = transcript_root / threshold / mode
            matrix_root = mode_root / "01_sample_selection/matrices"
            for specimen in SPECIMENS:
                if mode == "cell_separate":
                    specimen_root = matrix_root / specimen
                    for cell_root in sorted(path for path in specimen_root.iterdir() if path.is_dir()):
                        matches = sorted(cell_root.glob("*_after_cis_exclusion.csv"))
                        if len(matches) != 1:
                            raise ValueError(
                                f"Expected one checkpoint in {cell_root}, found {len(matches)}"
                            )
                        cell_type = cell_root.name
                        preprocess_root = mode_root / "03_preprocess" / specimen / cell_type
                        yield Network(
                            threshold,
                            mode,
                            specimen,
                            cell_type,
                            matches[0],
                            preprocess_root / "01_cut_diagnostic",
                            preprocess_root / "02_selected_cut",
                        )
                else:
                    checkpoint = (
                        matrix_root
                        / specimen
                        / f"TRANSCRIPTOMICS_{threshold}_{specimen}_CELL_INTEGRATED_after_cis_exclusion.csv"
                    )
                    preprocess_root = mode_root / "03_preprocess" / specimen
                    yield Network(
                        threshold,
                        mode,
                        specimen,
                        "cell_integrated",
                        checkpoint,
                        preprocess_root / "01_cut_diagnostic",
                        preprocess_root / "02_selected_cut",
                    )


def load_review(path: Path) -> dict[tuple[str, str, str, str], CutDecision]:
    review = pd.read_csv(path)
    required = {
        "threshold",
        "mode",
        "specimen",
        "cell_type",
        "decision",
        "confirmed",
        "selected_cut",
        "notes",
    }
    missing = required.difference(review.columns)
    if missing:
        raise ValueError(f"Review CSV is missing columns: {sorted(missing)}")

    decisions: dict[tuple[str, str, str, str], CutDecision] = {}
    for row in review.itertuples(index=False):
        key = (str(row.threshold), str(row.mode), str(row.specimen), str(row.cell_type))
        if key in decisions:
            raise ValueError(f"Duplicate review key: {key}")
        decision = str(row.decision).strip().lower()
        confirmed = str(row.confirmed).strip().lower() in {"yes", "true", "1"}
        if decision not in {"use", "reject"} or not confirmed:
            raise ValueError(f"Unresolved review row: {key}, decision={decision!r}")
        selected_cut: float | None = None
        if decision == "use":
            selected_cut = float(row.selected_cut)
            if not np.isfinite(selected_cut) or selected_cut <= 0:
                raise ValueError(f"Invalid selected cut for {key}: {row.selected_cut!r}")
        decisions[key] = CutDecision(
            decision=decision,
            selected_cut=selected_cut,
            source="review_csv",
            confirmed=True,
            notes="" if pd.isna(row.notes) else str(row.notes),
        )
    return decisions


def load_diagnostic_summary(network: Network) -> dict[str, Any]:
    path = network.diagnostic_root / "preprocess_summary.json"
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_decision(
    network: Network,
    reviewed: dict[tuple[str, str, str, str], CutDecision],
) -> tuple[CutDecision, dict[str, Any]]:
    summary = load_diagnostic_summary(network)
    if network.key in reviewed:
        return reviewed[network.key], summary
    flagged = int(summary.get("diagnostic_outlier_candidates", 0))
    if summary.get("status") == "complete_cut_diagnostic" and flagged != 0:
        raise ValueError(
            f"{network.key}: {flagged} candidates but no reviewed decision exists"
        )
    return CutDecision(
        decision="reject",
        selected_cut=None,
        source=(
            "automatic_no_cut_insufficient_samples"
            if summary.get("status") == "not_run_insufficient_samples"
            else "automatic_no_cut_zero_flags"
        ),
        confirmed=True,
        notes="No finite cut was reviewed because no diagnostic candidate was flagged.",
    ), summary


def plot_selected_dendrogram(
    linkage_matrix: np.ndarray,
    labels: list[str],
    decision: CutDecision,
    network_id: str,
    output_path: Path,
) -> list[int]:
    width = max(18.0, len(labels) / 10.0)
    figure, axis = plt.subplots(figsize=(width, 9), facecolor="white")
    color_threshold = decision.model_cut
    result = dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=6,
        color_threshold=color_threshold,
        above_threshold_color="#4C78A8",
        ax=axis,
    )
    if decision.cut_applied:
        axis.axhline(float(decision.selected_cut), color="#666666", linewidth=1.2, linestyle="--")
        cut_label = f"selected cut={float(decision.selected_cut):g}"
    else:
        cut_label = "no finite cut"
    axis.set_title(f"Transcriptomics sample clustering ({network_id}; {cut_label})")
    axis.set_xlabel("Samples")
    axis.set_ylabel("Euclidean distance; average linkage")
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return [int(value) for value in result["leaves"]]


def original_clusters(
    linkage_matrix: np.ndarray,
    decision: CutDecision,
) -> tuple[np.ndarray, int, dict[int, int]]:
    from PyWGCNA import WGCNA

    if not decision.cut_applied:
        labels = np.zeros(linkage_matrix.shape[0] + 1, dtype=int)
        return labels, 0, {0: int(labels.size)}
    labels = WGCNA.cutree(
        linkage_matrix, cutHeight=float(decision.selected_cut)
    )[:, 0].astype(int)
    values, counts = np.unique(labels, return_counts=True)
    cluster_counts = {int(label): int(count) for label, count in zip(values, counts)}
    largest_size = max(cluster_counts.values())
    largest_labels = [label for label, count in cluster_counts.items() if count == largest_size]
    if len(largest_labels) != 1:
        raise ValueError(
            f"Selected cut does not define a unique largest cluster: {cluster_counts}"
        )
    return labels, largest_labels[0], cluster_counts


@contextmanager
def pywgcna_retain_largest_cluster() -> Iterator[None]:
    """Make PyWGCNA's hard-coded label-0 rule retain the largest cluster."""
    from PyWGCNA import WGCNA

    original = WGCNA.cutree

    def relabelled_cutree(*args: Any, **kwargs: Any) -> np.ndarray:
        labels = np.asarray(original(*args, **kwargs)).copy()
        vector = labels[:, 0].astype(int)
        values, counts = np.unique(vector, return_counts=True)
        largest_size = int(counts.max())
        largest_labels = values[counts == largest_size]
        if largest_labels.size != 1:
            raise ValueError("PyWGCNA cut produced tied largest clusters; retention is ambiguous")
        largest = int(largest_labels[0])
        if largest != 0:
            swapped = vector.copy()
            swapped[vector == 0] = largest
            swapped[vector == largest] = 0
            labels[:, 0] = swapped
        return labels

    WGCNA.cutree = staticmethod(relabelled_cutree)
    try:
        yield
    finally:
        WGCNA.cutree = staticmethod(original)


def write_skipped_selected_cut(
    config: PipelineConfig,
    network: Network,
    decision: CutDecision,
    diagnostic_summary: dict[str, Any],
    review_path: Path,
) -> dict[str, Any]:
    network.selected_root.mkdir(parents=True, exist_ok=True)
    removed_features = network.selected_root / "removed_features.csv"
    removed_samples = network.selected_root / "removed_samples.csv"
    pd.DataFrame(columns=["gene_symbol", "removal_stage", "reason"]).to_csv(
        removed_features, index=False
    )
    pd.DataFrame(columns=["sample_id", "removal_stage", "reason"]).to_csv(
        removed_samples, index=False
    )
    summary_path = network.selected_root / "preprocess_summary.json"
    payload = {
        "status": "not_run_insufficient_samples",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "network_id": network.network_id,
        "reason": diagnostic_summary.get("reason", "diagnostic preprocess was skipped"),
        "decision": decision.decision,
        "decision_source": decision.source,
        "selected_cut": decision.serialized_cut,
        "findModules_run": False,
    }
    write_json(summary_path, payload)
    outputs = [removed_features, removed_samples, summary_path]
    write_run_manifest(
        network.selected_root / "run_manifest.json",
        pipeline_name="batch_corrected_transcriptomics_selected_cut_skipped",
        run_root=config.run_root,
        parameters={
            "network_id": network.network_id,
            "decision": decision.decision,
            "selected_cut": decision.serialized_cut,
        },
        generation_rules={
            "status": "explicitly skipped because default goodSamplesGenes requires at least 4 samples",
            "parameters_relaxed": False,
        },
        input_paths=[network.checkpoint_path, review_path],
        output_paths=outputs,
        code_paths=[Path(__file__)],
    )
    return payload


def run_network(
    config: PipelineConfig,
    network: Network,
    decision: CutDecision,
    diagnostic_summary: dict[str, Any],
    review_path: Path,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
) -> dict[str, Any]:
    from PyWGCNA import WGCNA

    if diagnostic_summary.get("status") == "not_run_insufficient_samples":
        return write_skipped_selected_cut(
            config, network, decision, diagnostic_summary, review_path
        )

    network.selected_root.mkdir(parents=True, exist_ok=True)
    input_matrix = read_expression_matrix(network.checkpoint_path)
    integrated_metadata = None
    if network.mode == "cell_integrated":
        _, integrated_metadata_path = integrated_paths(config, network.threshold)
        integrated_metadata = pd.read_csv(integrated_metadata_path)
    metadata = make_sample_metadata(
        input_matrix.index,
        clinical=clinical,
        cohort_map=cohort_map,
        fixed_cell_type=None if network.mode == "cell_integrated" else network.cell_type,
        integrated_metadata=integrated_metadata,
    )
    if metadata["diagnosis"].astype(str).eq("CIS").any():
        raise ValueError(f"{network.network_id}: CIS sample remains in checkpoint")
    sample_info = metadata.loc[
        input_matrix.index,
        ["project_id", "patient_id", "cohort", "specimen", "cell_type", "diagnosis", "sex", "age"],
    ].copy()

    model_name = f"TRANSCRIPTOMICS_{network.network_id}_selected_cut".replace(" ", "_")
    model = WGCNA(
        name=model_name,
        geneExp=input_matrix,
        sampleInfo=sample_info,
        species="Human",
        save=True,
        outputPath=network.selected_root.resolve().as_posix() + "/",
        figureType="png",
        cut=decision.model_cut,
    )

    tpm_keep = (input_matrix > model.TPMcutoff).any(axis=0).to_numpy()
    after_tpm = input_matrix.loc[:, tpm_keep]
    with pywgcna_pandas_compatibility():
        good_genes, good_samples, _ = WGCNA.goodSamplesGenes(after_tpm.T)
    good_genes = np.asarray(good_genes, dtype=bool)
    good_samples = np.asarray(good_samples, dtype=bool)
    before_cut = after_tpm.loc[good_samples, good_genes].copy()

    diagnostic_matrix_path = network.diagnostic_root / "matrix_before_cut.csv"
    diagnostic_matrix = read_expression_matrix(diagnostic_matrix_path)
    precut_equal = (
        list(before_cut.index) == list(diagnostic_matrix.index)
        and list(before_cut.columns) == list(diagnostic_matrix.columns)
        and np.array_equal(
            before_cut.to_numpy(dtype=float), diagnostic_matrix.to_numpy(dtype=float)
        )
    )
    if not precut_equal:
        raise ValueError(
            f"{network.network_id}: selected-cut pre-cut matrix differs from diagnostic"
        )

    linkage_matrix = WGCNA.hclust(
        pdist(before_cut.to_numpy(dtype=float), metric="euclidean"), method="average"
    )
    cluster_labels, kept_label, cluster_counts = original_clusters(linkage_matrix, decision)
    keep_mask = cluster_labels == kept_label
    label_zero_reassigned = bool(decision.cut_applied and kept_label != 0)

    with pywgcna_pandas_compatibility(), pywgcna_retain_largest_cluster():
        model.preprocess(show=False)
    actual = model.datExpr.to_df().copy()
    actual.index = actual.index.astype(str)
    actual.columns = actual.columns.astype(str)
    expected_after_cut = before_cut.iloc[keep_mask].copy()
    if (
        list(actual.index) != list(expected_after_cut.index)
        or list(actual.columns) != list(expected_after_cut.columns)
        or not np.array_equal(
            actual.to_numpy(dtype=float), expected_after_cut.to_numpy(dtype=float)
        )
    ):
        raise ValueError(f"{network.network_id}: selected-cut PyWGCNA output failed audit")

    model.analysis_provenance = {
        "stage": "02_selected_cut",
        "source_matrix": str(network.checkpoint_path),
        "diagnostic_matrix": str(diagnostic_matrix_path),
        "precut_matrix_equal_to_diagnostic": True,
        "cis_excluded": True,
        "transform": "none",
        "TPMcutoff": model.TPMcutoff,
        "goodSamplesGenes": (
            "PyWGCNA 2.2.1 default criteria with pandas compatibility fix; "
            "no threshold or inclusion rule changed"
        ),
        "decision": decision.decision,
        "decision_source": decision.source,
        "selected_cut": decision.serialized_cut,
        "cut_cluster_rule": "retain unique largest cluster",
        "original_kept_cluster_label": int(kept_label),
        "pywgcna_label_zero_reassigned": label_zero_reassigned,
        "dendrogram_distance": "Euclidean",
        "dendrogram_linkage": "average",
        "findModules_run": False,
    }
    model.saveWGCNA()
    object_path = network.selected_root / f"{model_name}.p"

    matrix_after_cut_path = network.selected_root / "matrix_after_cut.csv"
    write_expression_matrix(actual, matrix_after_cut_path, "sample_id")
    linkage_table = pd.DataFrame(
        linkage_matrix,
        columns=["left_node", "right_node", "distance", "cluster_size"],
    )
    linkage_table.insert(0, "merge_index", np.arange(len(linkage_table), dtype=int))
    linkage_table.to_csv(network.selected_root / "sample_linkage.csv", index=False)
    leaf_order = plot_selected_dendrogram(
        linkage_matrix,
        before_cut.index.astype(str).tolist(),
        decision,
        network.network_id,
        network.selected_root / "dendrogram_selected_cut.png",
    )
    pd.DataFrame({
        "leaf_position_1based": np.arange(1, len(leaf_order) + 1),
        "sample_id": before_cut.index[leaf_order].astype(str),
    }).to_csv(network.selected_root / "dendrogram_leaf_order.csv", index=False)

    assignments = sample_info.loc[before_cut.index].copy()
    assignments["original_cluster_label"] = cluster_labels
    assignments["kept_cluster_label"] = int(kept_label)
    assignments["selected_cut"] = (
        float(decision.selected_cut) if decision.cut_applied else np.nan
    )
    assignments["decision"] = decision.decision
    assignments["kept_by_selected_cut"] = keep_mask
    assignments["action"] = np.where(
        keep_mask, "retained", "removed_by_selected_cut"
    )
    assignments.reset_index(drop=True).to_csv(
        network.selected_root / "sample_cluster_assignments.csv", index=False
    )

    removed_sample_rows: list[dict[str, Any]] = []
    for sample_id in after_tpm.index[~good_samples]:
        removed_sample_rows.append({
            "sample_id": str(sample_id),
            "removal_stage": "goodSamplesGenes",
            "original_cluster_label": np.nan,
            "selected_cut": decision.serialized_cut,
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        })
    for sample_id, label in zip(
        before_cut.index[~keep_mask], cluster_labels[~keep_mask]
    ):
        removed_sample_rows.append({
            "sample_id": str(sample_id),
            "removal_stage": "dendrogram_selected_cut",
            "original_cluster_label": int(label),
            "selected_cut": decision.serialized_cut,
            "reason": "not in unique largest cluster at selected cut",
        })
    removed_samples = pd.DataFrame(
        removed_sample_rows,
        columns=[
            "sample_id",
            "removal_stage",
            "original_cluster_label",
            "selected_cut",
            "reason",
        ],
    )
    removed_samples.to_csv(network.selected_root / "removed_samples.csv", index=False)

    removed_feature_rows: list[dict[str, Any]] = []
    for feature in input_matrix.columns[~tpm_keep]:
        removed_feature_rows.append({
            "gene_symbol": str(feature),
            "removal_stage": "TPMcutoff",
            "reason": f"no sample value greater than TPMcutoff={model.TPMcutoff}",
        })
    for feature in after_tpm.columns[~good_genes]:
        removed_feature_rows.append({
            "gene_symbol": str(feature),
            "removal_stage": "goodSamplesGenes",
            "reason": "failed PyWGCNA 2.2.1 default criteria",
        })
    pd.DataFrame(
        removed_feature_rows,
        columns=["gene_symbol", "removal_stage", "reason"],
    ).to_csv(network.selected_root / "removed_features.csv", index=False)

    summary = {
        "status": "complete_selected_cut",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "network_id": network.network_id,
        "decision": decision.decision,
        "decision_source": decision.source,
        "decision_notes": decision.notes,
        "findModules_run": False,
        "input_samples_after_cis_exclusion": int(input_matrix.shape[0]),
        "input_genes": int(input_matrix.shape[1]),
        "TPMcutoff": float(model.TPMcutoff),
        "genes_removed_by_TPMcutoff": int((~tpm_keep).sum()),
        "genes_removed_by_goodSamplesGenes": int((~good_genes).sum()),
        "samples_removed_by_goodSamplesGenes": int((~good_samples).sum()),
        "selected_cut": decision.serialized_cut,
        "finite_cut_applied": decision.cut_applied,
        "cluster_counts": {str(key): value for key, value in cluster_counts.items()},
        "original_kept_cluster_label": int(kept_label),
        "kept_cluster_is_unique_largest": True,
        "pywgcna_label_zero_reassigned": label_zero_reassigned,
        "samples_removed_by_dendrogram_cut": int((~keep_mask).sum()),
        "samples_after_selected_cut": int(actual.shape[0]),
        "genes_after_selected_cut": int(actual.shape[1]),
        "precut_equal_to_cut_inf_diagnostic": True,
        "diagnostic_matrix_sha256": sha256(diagnostic_matrix_path),
        "matrix_after_cut_sha256": sha256(matrix_after_cut_path),
        "removed_sample_ids": removed_samples["sample_id"].astype(str).tolist(),
        "object": str(object_path),
    }
    summary_path = network.selected_root / "preprocess_summary.json"
    write_json(summary_path, summary)

    outputs = [
        matrix_after_cut_path,
        object_path,
        network.selected_root / "dendrogram_selected_cut.png",
        network.selected_root / "dendrogram_leaf_order.csv",
        network.selected_root / "sample_linkage.csv",
        network.selected_root / "sample_cluster_assignments.csv",
        network.selected_root / "removed_samples.csv",
        network.selected_root / "removed_features.csv",
        summary_path,
    ]
    write_run_manifest(
        network.selected_root / "run_manifest.json",
        pipeline_name="batch_corrected_transcriptomics_selected_cut",
        run_root=config.run_root,
        parameters={
            "network_id": network.network_id,
            "decision": decision.decision,
            "decision_source": decision.source,
            "selected_cut": decision.serialized_cut,
            "TPMcutoff": float(model.TPMcutoff),
            "goodSamplesGenes": "default criteria with pandas compatibility fix",
        },
        generation_rules={
            "execution": "new PyWGCNA object; full preprocess from CIS-excluded checkpoint",
            "input": "received batch-corrected values; no transformation or imputation",
            "precut_audit": "sample IDs, gene IDs, values, and hash equal cut=inf matrix_before_cut",
            "finite_cut_rule": "retain the unique largest cluster below the reviewed cut",
            "pywgcna_2_2_1_label_guard": (
                "swap largest cluster label with 0 inside preprocess only when scipy label 0 is not largest"
            ),
            "no_cut_rule": "reviewed reject and unflagged networks use cut=inf",
            "downstream": "findModules not run",
        },
        input_paths=[network.checkpoint_path, diagnostic_matrix_path, review_path],
        output_paths=outputs,
        code_paths=[
            Path(__file__),
            Path(__file__).parent / "run_transcriptomics_stage1.py",
            Path(__file__).parent / "pipeline_provenance.py",
        ],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )
    plt.close("all")
    del model
    gc.collect()
    return summary


def preflight(
    networks: list[Network],
    decisions: dict[tuple[str, str, str, str], tuple[CutDecision, dict[str, Any]]],
    reviewed: dict[tuple[str, str, str, str], CutDecision],
) -> pd.DataFrame:
    network_keys = {network.key for network in networks}
    unknown = set(reviewed).difference(network_keys)
    if unknown:
        raise ValueError(f"Review CSV contains unknown networks: {sorted(unknown)}")
    rows: list[dict[str, Any]] = []
    for network in networks:
        decision, diagnostic_summary = decisions[network.key]
        row: dict[str, Any] = {
            "threshold": network.threshold,
            "mode": network.mode,
            "specimen": network.specimen,
            "cell_type": network.cell_type,
            "decision": decision.decision,
            "decision_source": decision.source,
            "selected_cut": decision.selected_cut,
            "diagnostic_status": diagnostic_summary.get("status"),
            "diagnostic_flagged_count": int(
                diagnostic_summary.get("diagnostic_outlier_candidates", 0)
            ),
            "samples_before_cut": diagnostic_summary.get("matrix_before_cut_samples"),
            "original_label0_size": np.nan,
            "largest_cluster_size": np.nan,
            "original_kept_cluster_label": np.nan,
            "pywgcna_label_zero_reassigned": False,
        }
        if decision.cut_applied:
            linkage = pd.read_csv(network.diagnostic_root / "sample_linkage.csv")[
                ["left_node", "right_node", "distance", "cluster_size"]
            ].to_numpy(dtype=float)
            labels, kept_label, counts = original_clusters(linkage, decision)
            row.update({
                "original_label0_size": counts.get(0, 0),
                "largest_cluster_size": counts[kept_label],
                "original_kept_cluster_label": kept_label,
                "pywgcna_label_zero_reassigned": kept_label != 0,
                "samples_expected_after_cut": int((labels == kept_label).sum()),
                "samples_expected_removed_by_cut": int((labels != kept_label).sum()),
            })
        else:
            n = diagnostic_summary.get("matrix_before_cut_samples")
            row.update({
                "samples_expected_after_cut": n,
                "samples_expected_removed_by_cut": 0 if n is not None else np.nan,
            })
        rows.append(row)
    return pd.DataFrame(rows)


def run(config: PipelineConfig, *, preflight_only: bool = False) -> None:
    review_path = config.run_root / REVIEW_RELATIVE_PATH
    reviewed = load_review(review_path)
    networks = list(iter_networks(config))
    resolved = {
        network.key: resolve_decision(network, reviewed) for network in networks
    }
    audit = preflight(networks, resolved, reviewed)
    summary_root = config.experiment_root / "01_batch_corrected/transcriptomics"
    audit_path = summary_root / "selected_cut_preflight.csv"
    audit.to_csv(audit_path, index=False)
    if preflight_only:
        print(audit.to_string(index=False), flush=True)
        return

    clinical = pd.read_csv(config.clinical_path)
    cohort_map = build_source_cohort_map(config.run_root)
    execution_rows: list[dict[str, Any]] = []
    for position, network in enumerate(networks, start=1):
        decision, diagnostic_summary = resolved[network.key]
        print(
            f"[selected cut {position}/{len(networks)}] {network.network_id}: "
            f"decision={decision.decision}, cut={decision.serialized_cut}",
            flush=True,
        )
        summary = run_network(
            config,
            network,
            decision,
            diagnostic_summary,
            review_path,
            clinical,
            cohort_map,
        )
        execution_rows.append({
            "threshold": network.threshold,
            "mode": network.mode,
            "specimen": network.specimen,
            "cell_type": network.cell_type,
            "status": summary.get("status"),
            "decision": decision.decision,
            "decision_source": decision.source,
            "selected_cut": summary.get("selected_cut"),
            "samples_before_cut": summary.get("matrix_before_cut_samples"),
            "samples_removed_by_cut": summary.get("samples_removed_by_dendrogram_cut"),
            "samples_after_cut": summary.get("samples_after_selected_cut"),
            "genes_after_cut": summary.get("genes_after_selected_cut"),
            "pywgcna_label_zero_reassigned": summary.get(
                "pywgcna_label_zero_reassigned", False
            ),
        })

    execution_path = summary_root / "selected_cut_preprocess_summary.csv"
    pd.DataFrame(execution_rows).to_csv(execution_path, index=False)
    write_run_manifest(
        summary_root / "selected_cut_preprocess_manifest.json",
        pipeline_name="batch_corrected_transcriptomics_all_selected_cuts",
        run_root=config.run_root,
        parameters={
            "reviewed_use": sum(d.decision == "use" for d in reviewed.values()),
            "reviewed_reject": sum(d.decision == "reject" for d in reviewed.values()),
            "unreviewed_zero_flag_no_cut": sum(
                decision.source == "automatic_no_cut_zero_flags"
                for decision, _ in resolved.values()
            ),
            "insufficient_sample_networks": sum(
                summary.get("status") == "not_run_insufficient_samples"
                for _, summary in resolved.values()
            ),
        },
        generation_rules={
            "review_csv": "all rows must be confirmed use or reject",
            "use": "apply selected finite cut and retain unique largest cluster",
            "reject": "run preprocess with cut=inf",
            "absent_from_review": "allowed only when diagnostic flagged count is zero or preprocess was skipped",
            "stop_gate": "findModules not run",
        },
        input_paths=[review_path, audit_path],
        output_paths=[execution_path],
        code_paths=[Path(__file__)],
        packages=("numpy", "pandas", "matplotlib", "scipy", "PyWGCNA", "anndata"),
    )


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description="Apply reviewed transcriptomics cuts and stop before findModules."
    )
    parser.add_argument("--run-root", type=Path, default=default_root)
    parser.add_argument("--preflight-only", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = PipelineConfig(run_root=args.run_root.resolve())
    run(config, preflight_only=args.preflight_only)
    if not args.preflight_only:
        print("Completed reviewed transcriptomics selected-cut preprocess; findModules was not run.")


if __name__ == "__main__":
    main()
