"""Create direct disease-versus-HD module-eigengene association heatmaps.

PyWGCNA's multi-level categorical trait encoding compares each diagnosis with
all remaining diagnosis levels.  The scientific question for this project is
instead pairwise: keep only HD and one target disease, encode HD=0 and the
target disease=1, then correlate that binary indicator with every module
eigengene.  Pearson correlation with this binary indicator is exactly the
point-biserial correlation.

The existing official PyWGCNA one-hot outputs are intentionally preserved.
This script adds auditable pairwise-HD outputs beside them.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gc
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pointbiserialr

from pipeline_provenance import write_run_manifest
from run_official_kme_module_diagnosis import (
    CONDITION,
    PROTEOMICS_CONDITIONS,
    bh_fdr,
    find_model_used_by_findmodules,
    load_model,
    read_eigengenes,
    write_json,
)


REFERENCE_DIAGNOSIS = "HD"
TARGET_DIAGNOSES = ("MS", "MOGAD", "NMOSD")
MIN_GROUP_N_FOR_INFERENCE = 3
PAIRWISE_PREFIX = "module_diagnosis_pairwise_hd"


def finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def read_proper_modules(output_root: Path, eigengenes: pd.DataFrame) -> tuple[list[str], dict[str, int]]:
    """Return proper ME columns and module sizes from the audited kME summary."""

    summary_path = output_root / "module_kme_summary.csv"
    if not summary_path.exists():
        raise FileNotFoundError(
            f"Official kME summary is required before pairwise diagnosis analysis: {summary_path}"
        )
    summary = pd.read_csv(summary_path)
    required = {"module_color", "feature_count", "raw_unassigned_feature_count"}
    if not required.issubset(summary.columns):
        raise ValueError(f"Missing columns in {summary_path}: {sorted(required - set(summary.columns))}")

    fully_unassigned = summary["raw_unassigned_feature_count"].fillna(0).astype(int).eq(
        summary["feature_count"].astype(int)
    )
    proper_colors = summary.loc[~fully_unassigned, "module_color"].astype(str).tolist()
    size_map = dict(
        zip(
            summary["module_color"].astype(str),
            summary["feature_count"].astype(int),
            strict=True,
        )
    )
    proper_columns = [column for column in eigengenes.columns if column[2:] in proper_colors]
    missing = sorted(set("ME" + color for color in proper_colors) - set(proper_columns))
    if missing:
        raise ValueError(f"Proper module eigengenes are missing: {missing[:10]}")
    if not proper_columns:
        raise ValueError("No proper module eigengenes remain")
    return proper_columns, size_map


def calculate_pairwise_hd(
    eigengenes: pd.DataFrame,
    diagnosis: pd.Series,
    module_sizes: dict[str, int],
    *,
    network_id: str,
    omics: str,
    dimensions: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate target-disease-versus-HD point-biserial correlations."""

    diagnosis = diagnosis.astype("string").str.strip()
    diagnosis.index = diagnosis.index.astype(str)
    if set(eigengenes.index) != set(diagnosis.index):
        missing_metadata = sorted(set(eigengenes.index) - set(diagnosis.index))
        missing_mes = sorted(set(diagnosis.index) - set(eigengenes.index))
        raise ValueError(
            f"ME/metadata sample mismatch: missing_metadata={missing_metadata[:5]}, "
            f"missing_mes={missing_mes[:5]}"
        )
    diagnosis = diagnosis.reindex(eigengenes.index)

    contrasts = [f"{target} vs {REFERENCE_DIAGNOSIS}" for target in TARGET_DIAGNOSES]
    correlations = pd.DataFrame(np.nan, index=contrasts, columns=eigengenes.columns, dtype=float)
    pvalues = correlations.copy()
    qvalues = correlations.copy()
    statistics_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []

    for target in TARGET_DIAGNOSES:
        contrast = f"{target} vs {REFERENCE_DIAGNOSIS}"
        pair_mask = diagnosis.isin([REFERENCE_DIAGNOSIS, target])
        pair_diagnosis = diagnosis.loc[pair_mask]
        n_reference = int(pair_diagnosis.eq(REFERENCE_DIAGNOSIS).sum())
        n_target = int(pair_diagnosis.eq(target).sum())
        inference_ready = (
            n_reference >= MIN_GROUP_N_FOR_INFERENCE
            and n_target >= MIN_GROUP_N_FOR_INFERENCE
        )
        support_flag = (
            "INFERENCE_READY"
            if inference_ready
            else "NO_CASE"
            if min(n_reference, n_target) == 0
            else "SINGLE_CASE"
            if min(n_reference, n_target) == 1
            else "LOW_N"
        )
        support_rows.append(
            {
                "contrast": contrast,
                "reference_diagnosis": REFERENCE_DIAGNOSIS,
                "target_diagnosis": target,
                "n_reference": n_reference,
                "n_target": n_target,
                "n_pairwise": n_reference + n_target,
                "inference_ready": inference_ready,
                "support_flag": support_flag,
            }
        )

        binary = pair_diagnosis.eq(target).astype(float)
        for module in eigengenes.columns:
            color = module[2:] if module.startswith("ME") else module
            values = eigengenes.loc[pair_diagnosis.index, module].astype(float)
            valid = values.notna() & binary.notna()
            y = values.loc[valid]
            x = binary.loc[valid]
            valid_reference = int(x.eq(0).sum())
            valid_target = int(x.eq(1).sum())
            status = "OK"
            reason = ""
            correlation = np.nan
            p_value = np.nan
            if valid_reference == 0 or valid_target == 0:
                status = "MISSING_GROUP"
                reason = "HD or target diagnosis has no finite module eigengene value"
            elif y.nunique(dropna=True) < 2:
                status = "ZERO_VARIANCE_ME"
                reason = "Module eigengene is constant in the pairwise subset"
            else:
                result = pointbiserialr(x.to_numpy(dtype=float), y.to_numpy(dtype=float))
                correlation = float(result.statistic)
                p_value = float(result.pvalue)
                correlations.loc[contrast, module] = correlation
                pvalues.loc[contrast, module] = p_value

            hd_mean = finite_float(y.loc[x.eq(0)].mean())
            target_mean = finite_float(y.loc[x.eq(1)].mean())
            statistics_rows.append(
                {
                    "network_id": network_id,
                    "omics": omics,
                    **dimensions,
                    "contrast": contrast,
                    "reference_diagnosis": REFERENCE_DIAGNOSIS,
                    "target_diagnosis": target,
                    "module_eigengene": module,
                    "module_color": color,
                    "module_feature_count": module_sizes.get(color),
                    "n_reference": valid_reference,
                    "n_target": valid_target,
                    "n_pairwise": valid_reference + valid_target,
                    "inference_ready": (
                        valid_reference >= MIN_GROUP_N_FOR_INFERENCE
                        and valid_target >= MIN_GROUP_N_FOR_INFERENCE
                    ),
                    "support_flag": support_flag,
                    "point_biserial_r": correlation,
                    "p_value_two_sided": p_value,
                    "bh_fdr_q_within_network": np.nan,
                    "mean_me_hd": hd_mean,
                    "mean_me_target": target_mean,
                    "mean_difference_target_minus_hd": (
                        target_mean - hd_mean
                        if target_mean is not None and hd_mean is not None
                        else np.nan
                    ),
                    "calculation_status": status,
                    "non_estimable_reason": reason,
                }
            )

    statistics = pd.DataFrame(statistics_rows)
    eligible = statistics["inference_ready"].astype(bool) & statistics[
        "p_value_two_sided"
    ].notna()
    statistics.loc[eligible, "bh_fdr_q_within_network"] = bh_fdr(
        statistics.loc[eligible, "p_value_two_sided"].to_numpy(dtype=float)
    )
    for row in statistics.loc[eligible].itertuples(index=False):
        qvalues.loc[row.contrast, row.module_eigengene] = row.bh_fdr_q_within_network

    correlations.index.name = "contrast"
    pvalues.index.name = "contrast"
    qvalues.index.name = "contrast"
    return correlations, pvalues, qvalues, statistics, pd.DataFrame(support_rows)


def p_label(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def plot_pairwise_heatmap(
    correlations: pd.DataFrame,
    pvalues: pd.DataFrame,
    qvalues: pd.DataFrame,
    support: pd.DataFrame,
    module_sizes: dict[str, int],
    output_path: Path,
    title: str,
) -> None:
    """Plot signed pairwise effects with PyWGCNA-style r/p annotations."""

    annotations = pd.DataFrame("NA", index=correlations.index, columns=correlations.columns)
    support_lookup = support.set_index("contrast")
    for contrast in correlations.index:
        low_n = not bool(support_lookup.loc[contrast, "inference_ready"])
        for module in correlations.columns:
            r_value = correlations.loc[contrast, module]
            p_value = pvalues.loc[contrast, module]
            q_value = qvalues.loc[contrast, module]
            if np.isfinite(r_value):
                marker = "*" if np.isfinite(q_value) and q_value < 0.05 else ""
                low_n_marker = "†" if low_n else ""
                annotations.loc[contrast, module] = (
                    f"{r_value:+.2f}{marker}{low_n_marker}\n(p={p_label(p_value)})"
                )

    n_modules = correlations.shape[1]
    width = max(11.0, min(34.0, 0.72 * n_modules + 4.5))
    fig, ax = plt.subplots(figsize=(width, 5.8), constrained_layout=True)
    sns.heatmap(
        correlations,
        annot=annotations,
        fmt="",
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.6,
        linecolor="white",
        mask=correlations.isna(),
        cbar_kws={"label": "Point-biserial r (HD=0, target disease=1)", "shrink": 0.8},
        annot_kws={"fontsize": 7.5},
        ax=ax,
    )
    ax.set_facecolor("#E9EEF2")
    labels = []
    for module in correlations.columns:
        color = module[2:] if module.startswith("ME") else module
        labels.append(f"{color}\n(n={module_sizes.get(color, 'NA')})")
    ax.set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10, fontweight="bold")
    ax.set_xlabel("Final proper module (feature count)")
    ax.set_ylabel("Direct pairwise diagnosis contrast")
    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=14)
    ax.text(
        0,
        1.01,
        "Other diseases excluded per row · positive r = higher ME in target disease · "
        "* within-network BH q<0.05 · † group n<3",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="#405261",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def output_complete(output_root: Path) -> bool:
    required = (
        output_root / f"{PAIRWISE_PREFIX}_correlations.csv",
        output_root / f"{PAIRWISE_PREFIX}_pvalues.csv",
        output_root / f"{PAIRWISE_PREFIX}_qvalues.csv",
        output_root / f"{PAIRWISE_PREFIX}_statistics.csv",
        output_root / f"{PAIRWISE_PREFIX}_sample_counts.csv",
        output_root / "figures" / f"{PAIRWISE_PREFIX}_heatmap.png",
        output_root / "pairwise_hd_analysis_summary.json",
        output_root / "pairwise_hd_run_manifest.json",
    )
    return all(path.exists() and path.stat().st_size > 0 for path in required)


def run_network(
    *,
    run_root: Path,
    omics: str,
    network_id: str,
    source_model: Path,
    eigengenes_path: Path,
    output_root: Path,
    dimensions: dict[str, str],
    overwrite: bool,
) -> dict[str, Any]:
    if output_complete(output_root) and not overwrite:
        summary = pd.read_json(output_root / "pairwise_hd_analysis_summary.json", typ="series").to_dict()
        summary["resume_action"] = "reused_complete_output"
        return summary

    model = load_model(source_model)
    if "diagnosis" not in model.datExpr.obs.columns:
        raise ValueError(f"diagnosis is missing from {source_model}")
    eigengenes = read_eigengenes(eigengenes_path)
    eigengenes = eigengenes.reindex(model.datExpr.obs_names.astype(str))
    if eigengenes.isna().all(axis=1).any():
        raise ValueError("Module eigengenes do not cover all model samples")
    proper_columns, module_sizes = read_proper_modules(output_root, eigengenes)
    eigengenes = eigengenes[proper_columns]
    diagnosis = model.datExpr.obs["diagnosis"].copy()
    diagnosis.index = model.datExpr.obs_names.astype(str)

    correlations, pvalues, qvalues, statistics, support = calculate_pairwise_hd(
        eigengenes,
        diagnosis,
        module_sizes,
        network_id=network_id,
        omics=omics,
        dimensions=dimensions,
    )
    paths = {
        "correlations": output_root / f"{PAIRWISE_PREFIX}_correlations.csv",
        "pvalues": output_root / f"{PAIRWISE_PREFIX}_pvalues.csv",
        "qvalues": output_root / f"{PAIRWISE_PREFIX}_qvalues.csv",
        "statistics": output_root / f"{PAIRWISE_PREFIX}_statistics.csv",
        "support": output_root / f"{PAIRWISE_PREFIX}_sample_counts.csv",
        "heatmap": output_root / "figures" / f"{PAIRWISE_PREFIX}_heatmap.png",
    }
    correlations.to_csv(paths["correlations"])
    pvalues.to_csv(paths["pvalues"])
    qvalues.to_csv(paths["qvalues"])
    statistics.to_csv(paths["statistics"], index=False)
    support.to_csv(paths["support"], index=False)

    display_bits = [omics.capitalize(), dimensions.get("specimen", "ALL")]
    if dimensions.get("cell_type") and dimensions["cell_type"] not in ("", "nan", "ALL"):
        display_bits.append(dimensions["cell_type"])
    if dimensions.get("strategy"):
        display_bits.append(dimensions["strategy"])
    plot_pairwise_heatmap(
        correlations,
        pvalues,
        qvalues,
        support,
        module_sizes,
        paths["heatmap"],
        " · ".join(display_bits) + " — module–diagnosis association",
    )

    finite_p = statistics["p_value_two_sided"].dropna()
    finite_q = statistics["bh_fdr_q_within_network"].dropna()
    summary = {
        "status": "complete",
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
        "omics": omics,
        "network_id": network_id,
        **dimensions,
        "samples": int(model.datExpr.n_obs),
        "proper_modules": int(eigengenes.shape[1]),
        "diagnosis_levels_observed": sorted(diagnosis.dropna().astype(str).unique().tolist()),
        "contrasts": support["contrast"].tolist(),
        "tests_total": int(len(statistics)),
        "tests_estimable": int(statistics["p_value_two_sided"].notna().sum()),
        "tests_inference_ready": int(statistics["inference_ready"].astype(bool).sum()),
        "minimum_p_value": float(finite_p.min()) if len(finite_p) else None,
        "minimum_bh_fdr_q_within_network": float(finite_q.min()) if len(finite_q) else None,
        "significant_tests_q_lt_0_05": int(finite_q.lt(0.05).sum()),
        "method": "pairwise point-biserial correlation; HD=0, target disease=1; other diseases excluded",
        "heatmap_annotations": "signed r and two-sided unadjusted p; star marks within-network BH q<0.05",
        "inference_ready_rule": (
            f"n_HD >= {MIN_GROUP_N_FOR_INFERENCE} and n_target >= {MIN_GROUP_N_FOR_INFERENCE}"
        ),
        "non_estimable_policy": "missing group or zero-variance ME remains NA; never replace with zero",
    }
    summary_path = output_root / "pairwise_hd_analysis_summary.json"
    write_json(summary_path, summary)
    manifest_path = output_root / "pairwise_hd_run_manifest.json"
    write_run_manifest(
        manifest_path,
        pipeline_name=f"{omics}_pairwise_hd_module_diagnosis",
        run_root=run_root,
        parameters={
            "reference_diagnosis": REFERENCE_DIAGNOSIS,
            "target_diagnoses": list(TARGET_DIAGNOSES),
            "minimum_group_n_for_inference": MIN_GROUP_N_FOR_INFERENCE,
            "alternative": "two-sided",
            **dimensions,
        },
        generation_rules={
            "sample_subset": "for each row keep only HD and the target diagnosis; exclude every other diagnosis",
            "diagnosis_encoding": "HD=0; target diagnosis=1",
            "statistic": "scipy.stats.pointbiserialr; mathematically Pearson r with a binary diagnosis indicator",
            "direction": "positive r means higher module eigengene in target diagnosis than HD",
            "multiple_testing": "BH-FDR across inference-ready module x target-diagnosis tests within this network",
            "display": "signed r with unadjusted two-sided p; star marks within-network q<0.05; low-support estimates carry dagger",
            "missingness": "missing group or zero-variance eigengene is shown as NA, never zero",
            "module_scope": "saved final proper module eigengenes; fully raw-unassigned eigengene excluded",
        },
        input_paths=[source_model, eigengenes_path, output_root / "module_kme_summary.csv"],
        output_paths=[*paths.values(), summary_path],
        code_paths=[Path(__file__), Path(__file__).parent / "run_official_kme_module_diagnosis.py"],
        packages=("numpy", "pandas", "matplotlib", "scipy", "seaborn", "PyWGCNA", "anndata"),
    )
    del model, eigengenes, correlations, pvalues, qvalues, statistics
    gc.collect()
    return summary


def run_proteomics(run_root: Path, overwrite: bool, limit: int | None) -> list[dict[str, Any]]:
    experiment = run_root / "experiments/01_batch_corrected/proteomics"
    rows: list[dict[str, Any]] = []
    for specimen in ("CSF", "PB"):
        source_model = experiment / "03_preprocess" / specimen / "02_selected_cut" / f"PROT_{specimen}_selected_cut.p"
        for condition in PROTEOMICS_CONDITIONS:
            if limit is not None and len(rows) >= limit:
                return rows
            module_root = experiment / "05_find_modules" / specimen / condition
            output_root = experiment / "06_kme_module_diagnosis" / specimen / condition
            print(f"[proteomics pairwise-HD] {specimen}/{condition}", flush=True)
            rows.append(
                run_network(
                    run_root=run_root,
                    omics="proteomics",
                    network_id=f"proteomics_{specimen}_{condition}",
                    source_model=source_model,
                    eigengenes_path=module_root / "module_eigengenes.csv",
                    output_root=output_root,
                    dimensions={"specimen": specimen, "condition": condition},
                    overwrite=overwrite,
                )
            )
    return rows


def run_transcriptomics(run_root: Path, overwrite: bool, limit: int | None) -> list[dict[str, Any]]:
    experiment = run_root / "experiments/04_variance_filtered/transcriptomics"
    registry = pd.read_csv(experiment / "find_modules_summary.csv")
    registry = registry.loc[registry["status"].eq("complete_findModules")].copy()
    rows: list[dict[str, Any]] = []
    for row in registry.itertuples(index=False):
        if limit is not None and len(rows) >= limit:
            break
        network_root = experiment / row.threshold / row.strategy / row.mode / row.specimen
        if row.mode == "cell_separate":
            network_root = network_root / row.cell_type
        module_root = network_root / "05_find_modules" / CONDITION
        output_root = network_root / "06_kme_module_diagnosis" / CONDITION
        source_model = find_model_used_by_findmodules(
            run_root, module_root, network_root / "03_preprocess" / "02_selected_cut"
        )
        print(f"[transcriptomics pairwise-HD] {row.network_id} / {row.strategy}", flush=True)
        rows.append(
            run_network(
                run_root=run_root,
                omics="transcriptomics",
                network_id=str(row.network_id),
                source_model=source_model,
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
    return rows


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--omics", choices=("all", "proteomics", "transcriptomics"), default="all")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_root = args.run_root.resolve()
    if args.omics in ("all", "proteomics"):
        rows = run_proteomics(run_root, args.overwrite, args.limit)
        write_summary(
            run_root / "experiments/01_batch_corrected/proteomics/06_kme_module_diagnosis/pairwise_hd_summary.csv",
            rows,
        )
    if args.omics in ("all", "transcriptomics"):
        rows = run_transcriptomics(run_root, args.overwrite, args.limit)
        write_summary(
            run_root / "experiments/04_variance_filtered/transcriptomics/pairwise_hd_module_diagnosis_summary.csv",
            rows,
        )


if __name__ == "__main__":
    main()
