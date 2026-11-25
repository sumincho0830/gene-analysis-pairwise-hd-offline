"""Shared transcriptomics QC utilities for the results_20261119 pipeline.

All calculations use the stored expression values without log transformation,
feature scaling, missing-value imputation, or top-N feature selection. PCA
removes an entire gene if any value is non-finite in the analysis scope and
then removes zero-variance genes; the retained values are mean-centred only by
``sklearn.decomposition.PCA``.
"""

from __future__ import annotations

from dataclasses import dataclass
import html
import json
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


PROJECT_RE = re.compile(r"^(PRJ-[A-Za-z]+-\d+)")
SPECIMEN_RE = re.compile(r"_(CSF|PB)(?:\d+)?(?:_|$)")
CELL_SUFFIX_RE = re.compile(r"__CELL__(?:__)?(?P<cell>.*)$")
OUTLIER_Z_THRESHOLD = 5.0
SPECIMEN_ORDER = {"CSF": 0, "PB": 1, "Unknown": 2}
MARKERS = ("o", "s", "^", "D", "P", "X", "v", "<", ">", "h")
@dataclass(frozen=True)
class QCOutputs:
    sample_statistics: pd.DataFrame
    gene_statistics: pd.DataFrame | None
    pca_scores: pd.DataFrame
    pca_features_used: int
    paths: tuple[Path, ...]


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def read_expression_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if frame.shape[1] < 2:
        raise ValueError(f"Expression matrix has no feature columns: {path}")
    sample_ids = frame.iloc[:, 0].astype(str)
    if sample_ids.duplicated().any():
        duplicated = sample_ids[sample_ids.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate sample IDs in {path}: {duplicated}")
    numeric = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    invalid = frame.iloc[:, 1:].notna() & numeric.isna()
    if invalid.to_numpy().any():
        raise ValueError(f"Non-numeric expression values in {path}")
    values = numeric.to_numpy(dtype=float, copy=False)
    if np.isinf(values).any():
        raise ValueError(f"Infinite expression values in {path}")
    numeric.index = pd.Index(sample_ids, name="sample_id")
    numeric.columns = numeric.columns.astype(str)
    return numeric.astype(float)


def write_expression_matrix(matrix: pd.DataFrame, path: Path, index_label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, index=True, index_label=index_label)


def extract_specimen(sample_id: str) -> str:
    clean = str(sample_id).split("__CELL__", 1)[0]
    match = SPECIMEN_RE.search(clean)
    return match.group(1) if match else "Unknown"


def build_source_cohort_map(run_root: Path) -> dict[str, str]:
    root = run_root / "data/transcriptomics/cohorts"
    mapping: dict[str, str] = {}
    for path in sorted(root.glob("*/*.csv")):
        sample_id = path.stem
        cohort = path.parent.name
        previous = mapping.get(sample_id)
        if previous is not None and previous != cohort:
            raise ValueError(f"Sample appears in two transcriptomics cohorts: {sample_id}")
        mapping[sample_id] = cohort
    if not mapping:
        raise FileNotFoundError(f"No source transcriptomics files under {root}")
    return mapping


def make_sample_metadata(
    sample_ids: pd.Index,
    *,
    clinical: pd.DataFrame,
    cohort_map: dict[str, str],
    fixed_cell_type: str | None = None,
    integrated_metadata: pd.DataFrame | None = None,
) -> pd.DataFrame:
    integrated_lookup: pd.DataFrame | None = None
    if integrated_metadata is not None:
        key = integrated_metadata.columns[0]
        integrated_lookup = integrated_metadata.copy()
        integrated_lookup[key] = integrated_lookup[key].astype(str)
        integrated_lookup = integrated_lookup.set_index(key, drop=False)

    records: list[dict[str, Any]] = []
    for row_id in sample_ids.astype(str):
        if integrated_lookup is not None:
            if row_id not in integrated_lookup.index:
                raise ValueError(f"Integrated metadata is missing row {row_id}")
            integrated_row = integrated_lookup.loc[row_id]
            source_sample_id = str(integrated_row["source_sample_id"])
            cell_type = str(integrated_row["cell_type"])
        else:
            source_sample_id = row_id.split("__CELL__", 1)[0]
            cell_type = fixed_cell_type or "Unknown"

        project_match = PROJECT_RE.match(source_sample_id)
        if project_match is None:
            raise ValueError(f"Cannot parse project ID from {source_sample_id}")
        project_id = project_match.group(1)
        cohort = cohort_map.get(source_sample_id)
        if cohort is None:
            raise ValueError(f"No source cohort mapping for {source_sample_id}")
        candidates = clinical.loc[
            clinical["Project_ID"].astype(str).eq(project_id)
            & clinical["Cohort"].astype(str).eq(cohort)
        ]
        if len(candidates) != 1:
            raise ValueError(
                f"Expected one clinical row for ({project_id}, {cohort}), found {len(candidates)}"
            )
        clinical_row = candidates.iloc[0]
        records.append({
            "sample_id": row_id,
            "source_sample_id": source_sample_id,
            "project_id": project_id,
            "patient_id": clinical_row["Patient_ID"],
            "cohort": cohort,
            "specimen": extract_specimen(source_sample_id),
            "cell_type": cell_type,
            "diagnosis": clinical_row["Diagnosis"],
            "sex": clinical_row["Sex"],
            "age": clinical_row["Age"],
        })
    metadata = pd.DataFrame(records).set_index("sample_id", drop=False)
    if metadata.index.duplicated().any():
        raise ValueError("Duplicate metadata row IDs")
    return metadata.loc[sample_ids.astype(str)]


def robust_z(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    median = numeric.median()
    mad = (numeric - median).abs().median()
    result = pd.Series(np.nan, index=numeric.index, dtype=float)
    if np.isfinite(mad) and mad > 0:
        result = 0.6744897501960817 * (numeric - median) / mad
    return result


def calculate_sample_statistics(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    finite = np.isfinite(values)
    finite_values = np.where(finite, values, np.nan)
    q25 = np.nanpercentile(finite_values, 25, axis=1)
    q75 = np.nanpercentile(finite_values, 75, axis=1)
    table = metadata.copy()
    table["feature_count"] = matrix.shape[1]
    table["observed_feature_count"] = finite.sum(axis=1)
    table["missing_feature_count"] = (~finite).sum(axis=1)
    table["missing_fraction"] = (~finite).mean(axis=1)
    table["zero_fraction"] = np.divide(
        ((values == 0) & finite).sum(axis=1),
        finite.sum(axis=1),
        out=np.full(values.shape[0], np.nan, dtype=float),
        where=finite.sum(axis=1) > 0,
    )
    table["sample_mean"] = np.nanmean(finite_values, axis=1)
    table["sample_median"] = np.nanmedian(finite_values, axis=1)
    table["sample_std"] = np.nanstd(finite_values, axis=1, ddof=0)
    table["sample_iqr"] = q75 - q25
    table["sample_median_robust_z"] = robust_z(table["sample_median"])
    table["sample_iqr_robust_z"] = robust_z(table["sample_iqr"])
    table["zero_fraction_robust_z"] = robust_z(table["zero_fraction"])
    descriptive_columns = [
        "sample_median_robust_z",
        "sample_iqr_robust_z",
        "zero_fraction_robust_z",
    ]
    table["descriptive_outlier_candidate"] = table[descriptive_columns].abs().ge(
        OUTLIER_Z_THRESHOLD
    ).any(axis=1)
    return table


def calculate_gene_statistics(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(dtype=float)
    finite = np.isfinite(values)
    finite_values = np.where(finite, values, np.nan)
    q25 = np.nanpercentile(finite_values, 25, axis=0)
    q75 = np.nanpercentile(finite_values, 75, axis=0)
    return pd.DataFrame({
        "gene_symbol": matrix.columns.astype(str),
        "sample_count": matrix.shape[0],
        "observed_sample_count": finite.sum(axis=0),
        "missing_sample_count": (~finite).sum(axis=0),
        "missing_fraction": (~finite).mean(axis=0),
        "zero_fraction": np.divide(
            ((values == 0) & finite).sum(axis=0),
            finite.sum(axis=0),
            out=np.full(values.shape[1], np.nan, dtype=float),
            where=finite.sum(axis=0) > 0,
        ),
        "mean": np.nanmean(finite_values, axis=0),
        "median": np.nanmedian(finite_values, axis=0),
        "std": np.nanstd(finite_values, axis=0, ddof=0),
        "iqr": q75 - q25,
        "minimum": np.nanmin(finite_values, axis=0),
        "maximum": np.nanmax(finite_values, axis=0),
    })


def calculate_pca(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, int, int]:
    values = matrix.to_numpy(dtype=float)
    complete = np.isfinite(values).all(axis=0)
    complete_values = values[:, complete]
    nonconstant = np.nanvar(complete_values, axis=0, ddof=0) > 0
    pca_values = complete_values[:, nonconstant]
    if pca_values.shape[1] == 0:
        raise ValueError("No complete-case nonconstant genes remain for PCA")
    components = min(3, pca_values.shape[0] - 1, pca_values.shape[1])
    if components < 2:
        raise ValueError("At least three samples are required for PCA figures")
    solver = "randomized" if min(pca_values.shape) > 50 else "full"
    pca = PCA(n_components=components, svd_solver=solver, random_state=0)
    scores = pca.fit_transform(pca_values)
    table = metadata.copy()
    for index in range(3):
        table[f"PC{index + 1}"] = scores[:, index] if index < components else 0.0
        table[f"PC{index + 1}_explained_variance_ratio"] = (
            float(pca.explained_variance_ratio_[index]) if index < components else 0.0
        )
    distance = np.sqrt(np.sum(scores[:, :components] ** 2, axis=1))
    table["pca_distance_pc1_pc3"] = distance
    table["pca_distance_robust_z"] = robust_z(pd.Series(distance, index=table.index))
    table["pca_outlier_candidate"] = table["pca_distance_robust_z"].ge(
        OUTLIER_Z_THRESHOLD
    )
    removed_missing = int((~complete).sum())
    removed_constant = int((~nonconstant).sum())
    return table, int(pca_values.shape[1]), removed_missing + removed_constant


def ordered_ids(metadata: pd.DataFrame) -> list[str]:
    order = metadata.copy()
    order["_specimen_rank"] = order["specimen"].map(SPECIMEN_ORDER).fillna(99)
    order["_sample_id_sort"] = order.index.astype(str)
    return order.sort_values(
        ["_specimen_rank", "cohort", "cell_type", "source_sample_id", "_sample_id_sort"]
    ).index.astype(str).tolist()


def _category_colors(values: pd.Series) -> dict[str, str]:
    categories = sorted(values.fillna("Unknown").astype(str).unique())
    cmap = plt.get_cmap("tab20")
    return {value: mcolors.to_hex(cmap(index % 20)) for index, value in enumerate(categories)}


def _outlier_mask(frame: pd.DataFrame) -> pd.Series:
    """Return the saved PCA candidate flag without recalculating it."""
    column = "pca_outlier_candidate" if "pca_outlier_candidate" in frame else "pca_outlier"
    return frame[column].fillna(False).astype(bool)


def _ranked_outliers(frame: pd.DataFrame) -> pd.DataFrame:
    """Order saved candidates only for readable text placement in figures."""
    outliers = frame.loc[_outlier_mask(frame)].copy()
    sort_columns = [
        column
        for column in (
            "pca_distance_robust_z",
            "pca_robust_distance_pc1_pc5",
            "pca_distance_pc1_pc3",
        )
        if column in outliers
    ]
    if sort_columns:
        outliers = outliers.sort_values(sort_columns, ascending=False, na_position="last")
    if "source_sample_id" in outliers:
        outliers["_outlier_label"] = outliers["source_sample_id"].astype(str)
    elif "sample_id" in outliers:
        outliers["_outlier_label"] = outliers["sample_id"].astype(str)
    else:
        outliers["_outlier_label"] = outliers.index.astype(str)
    duplicated = outliers["_outlier_label"].duplicated(keep=False)
    if duplicated.any() and "cell_type" in outliers:
        outliers.loc[duplicated, "_outlier_label"] = (
            outliers.loc[duplicated, "_outlier_label"]
            + " ["
            + outliers.loc[duplicated, "cell_type"].fillna("Unknown").astype(str)
            + "]"
        )
    return outliers


def plot_intensity_distribution(
    matrix: pd.DataFrame,
    sample_statistics: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    values = matrix.to_numpy(dtype=float).ravel()
    values = values[np.isfinite(values)]
    if len(values) > 1_000_000:
        rng = np.random.default_rng(0)
        values = rng.choice(values, size=1_000_000, replace=False)
    figure, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor="white")
    axes[0].hist(values, bins=100, color="#4C78A8", alpha=0.85)
    axes[0].set_xlabel("Stored expression value")
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Expression-value distribution")
    groups = []
    labels = []
    for specimen in ("CSF", "PB"):
        group = sample_statistics.loc[
            sample_statistics["specimen"].eq(specimen), "sample_median"
        ].dropna()
        if len(group):
            groups.append(group.to_numpy())
            labels.append(specimen)
    if not groups:
        groups = [sample_statistics["sample_median"].dropna().to_numpy()]
        labels = ["ALL"]
    axes[1].boxplot(groups, tick_labels=labels, showfliers=True)
    axes[1].set_ylabel("Sample median stored expression")
    axes[1].set_title("Sample-level intensity")
    figure.suptitle(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_heatmap(
    matrix: pd.DataFrame,
    sample_statistics: pd.DataFrame,
    output_path: Path,
    title: str,
) -> None:
    order = ordered_ids(sample_statistics)
    visual = matrix.loc[order].to_numpy(dtype=float)
    missing = ~np.isfinite(visual)
    finite_values = visual[np.isfinite(visual)]
    if len(finite_values):
        vmin, vmax = np.percentile(finite_values, [1, 99])
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
            vmin, vmax = float(np.nanmin(finite_values)), float(np.nanmax(finite_values) + 1e-12)
    else:
        vmin, vmax = 0.0, 1.0
    height = max(8.0, min(18.0, len(order) / 20 + 6.0))
    figure = plt.figure(figsize=(24, height), facecolor="white")
    grid = figure.add_gridspec(2, 2, width_ratios=[30, 2], height_ratios=[1, 4], hspace=0.18, wspace=0.03)
    ax_missing = figure.add_subplot(grid[0, 0])
    ax_value = figure.add_subplot(grid[1, 0])
    ax_flag = figure.add_subplot(grid[:, 1])
    missing_cmap = mcolors.ListedColormap(["#24384A", "#F4A261"])
    ax_missing.imshow(missing.astype(int), aspect="auto", interpolation="nearest", cmap=missing_cmap, vmin=0, vmax=1)
    ax_missing.set_title(f"Missingness ({missing.mean() * 100:.2f}% of cells)")
    ax_missing.set_xticks([])
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("lightgrey")
    image = ax_value.imshow(np.ma.masked_invalid(visual), aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
    ax_value.set_xlabel("Genes (stored order)")
    ax_value.set_ylabel("Samples")
    ax_value.set_title("Stored expression intensity (within-matrix 1st–99th percentile)")
    colorbar = figure.colorbar(image, ax=ax_value, fraction=0.015, pad=0.01)
    colorbar.set_label("Stored expression value")
    ordered_stats = sample_statistics.loc[order]
    label_columns = ["specimen", "cohort"]
    if ordered_stats["cell_type"].astype(str).nunique() > 1:
        label_columns.append("cell_type")
    group_labels = ordered_stats[label_columns].astype(str).agg(" / ".join, axis=1)
    group_boundaries = np.flatnonzero(
        group_labels.to_numpy()[1:] != group_labels.to_numpy()[:-1]
    ) + 0.5
    specimen_values = ordered_stats["specimen"].astype(str).to_numpy()
    specimen_boundaries = np.flatnonzero(specimen_values[1:] != specimen_values[:-1]) + 0.5
    group_ticks: list[float] = []
    group_tick_labels: list[str] = []
    start = 0
    for index in range(1, len(group_labels) + 1):
        if index == len(group_labels) or group_labels.iloc[index] != group_labels.iloc[start]:
            group_ticks.append((start + index - 1) / 2)
            group_tick_labels.append(group_labels.iloc[start])
            start = index
    for axis in (ax_missing, ax_value):
        for boundary in group_boundaries:
            axis.axhline(boundary, color="white", linewidth=0.7)
        for boundary in specimen_boundaries:
            axis.axhline(boundary, color="#00A6A6", linewidth=1.2)
        axis.set_yticks(group_ticks, group_tick_labels, fontsize=6)
    flagged = ordered_stats["outlier_candidate"].fillna(False).to_numpy(dtype=bool)
    flag_image = np.where(flagged[:, None], 1, 0)
    ax_flag.imshow(flag_image, aspect="auto", interpolation="nearest", cmap=mcolors.ListedColormap(["white", "red"]), vmin=0, vmax=1)
    ax_flag.set_xticks([])
    positions = np.flatnonzero(flagged)
    ax_flag.set_yticks(positions)
    ax_flag.set_yticklabels(ordered_stats.iloc[positions]["source_sample_id"], fontsize=6, color="red")
    ax_flag.yaxis.tick_right()
    ax_flag.set_title("Outlier", fontsize=9)
    figure.suptitle(title, y=0.995)
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _scatter_2d(
    axis: plt.Axes,
    frame: pd.DataFrame,
    color_column: str,
    marker_column: str,
) -> None:
    color_values = frame[color_column].fillna("Unknown").astype(str)
    marker_values = frame[marker_column].fillna("Unknown").astype(str)
    colors = _category_colors(color_values)
    markers = {
        value: MARKERS[index % len(MARKERS)]
        for index, value in enumerate(sorted(marker_values.unique()))
    }
    plotting = frame.assign(_color_value=color_values, _marker_value=marker_values)
    for (color_value, marker_value), group in plotting.groupby(
        ["_color_value", "_marker_value"], dropna=False
    ):
        axis.scatter(
            group["PC1"], group["PC2"],
            c=colors[color_value], marker=markers[marker_value],
            s=38, alpha=0.85, edgecolors="black", linewidths=0.25,
        )
    outliers = _ranked_outliers(frame)
    if len(outliers):
        axis.scatter(
            outliers["PC1"], outliers["PC2"], facecolors="none",
            edgecolors="#E15759", s=120, linewidths=1.6,
        )
        offsets = ((5, 5), (5, 13), (5, -11), (-5, 5), (-5, 13), (-5, -11))
        for rank, (_, row) in enumerate(outliers.iterrows()):
            dx, dy = offsets[rank % len(offsets)]
            axis.annotate(
                row["_outlier_label"], (row["PC1"], row["PC2"]),
                xytext=(dx, dy), textcoords="offset points", fontsize=6,
                ha="left" if dx > 0 else "right",
                arrowprops={"arrowstyle": "-", "color": "#777777", "lw": 0.35},
            )
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor=colors[value],
            markeredgecolor="none", label=value, markersize=6,
        )
        for value in sorted(color_values.unique())
    ]
    handles.extend(
        Line2D(
            [0], [0], marker=markers[value], linestyle="", color="#444444",
            label=value, markersize=6,
        )
        for value in sorted(marker_values.unique())
    )
    if len(outliers):
        handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="", markerfacecolor="none",
                markeredgecolor="#E15759", label="PCA outlier candidate", markersize=8,
            )
        )
    axis.set_title(f"Color: {color_column}; marker: {marker_column}")
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    axis.grid(alpha=0.2)
    axis.legend(handles=handles, loc="best", frameon=False, fontsize=6, ncol=2)


def plot_pca_metadata(frame: pd.DataFrame, output_path: Path, title: str) -> None:
    views = [
        ("cohort", "diagnosis"),
        ("diagnosis", "cohort"),
        ("cell_type", "specimen"),
        ("specimen", "cell_type"),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(16, 13), facecolor="white")
    for axis, (color_column, marker_column) in zip(axes.ravel(), views):
        _scatter_2d(axis, frame, color_column, marker_column)
    figure.suptitle(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_pca_3d(frame: pd.DataFrame, output_path: Path, title: str) -> None:
    figure = plt.figure(figsize=(12, 10), facecolor="white")
    axis = figure.add_subplot(111, projection="3d")
    colors = _category_colors(frame["cohort"])
    for cohort, group in frame.groupby("cohort", dropna=False):
        axis.scatter(group["PC1"], group["PC2"], group["PC3"], c=colors[str(cohort)], label=str(cohort), s=34, alpha=0.8)
    outliers = _ranked_outliers(frame)
    if len(outliers):
        axis.scatter(
            outliers["PC1"], outliers["PC2"], outliers["PC3"],
            facecolors="none", edgecolors="#E15759", s=130, linewidths=1.5,
            label="PCA outlier candidate", depthshade=False,
        )
        for _, row in outliers.iterrows():
            axis.text(
                row["PC1"], row["PC2"], row["PC3"],
                row["_outlier_label"], fontsize=6,
            )
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    axis.set_zlabel("PC3")
    axis.set_title(title)
    axis.legend(title="Cohort", fontsize=7, loc="best")
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_pca_3d_interactive(
    frame: pd.DataFrame,
    output_path: Path,
    title: str,
    *,
    include_plotlyjs: bool | str = True,
) -> None:
    import plotly.express as px

    hover = [
        column for column in (
            "source_sample_id", "specimen", "cell_type", "cohort", "diagnosis",
            "sex", "age", "pca_distance_robust_z", "pca_outlier_candidate",
        ) if column in frame.columns
    ]
    figure = px.scatter_3d(
        frame.reset_index(drop=True),
        x="PC1", y="PC2", z="PC3",
        color="cohort", symbol="diagnosis",
        hover_name="source_sample_id", hover_data=hover,
        title=title,
    )
    figure.update_traces(marker={"size": 5, "opacity": 0.8})
    outliers = _ranked_outliers(frame)
    if len(outliers):
        import plotly.graph_objects as go

        text = outliers["_outlier_label"].tolist()
        hover_text = [
            "<br>".join(
                [
                    f"<b>{row['_outlier_label']}</b>",
                    f"Specimen: {row.get('specimen', 'Unknown')}",
                    f"Cohort: {row.get('cohort', 'Unknown')}",
                    f"Cell type: {row.get('cell_type', 'Unknown')}",
                    f"Diagnosis: {row.get('diagnosis', 'Unknown')}",
                    f"PCA distance robust z: {row.get('pca_distance_robust_z', np.nan):.4g}",
                ]
            )
            for _, row in outliers.iterrows()
        ]
        figure.add_trace(
            go.Scatter3d(
                x=outliers["PC1"], y=outliers["PC2"], z=outliers["PC3"],
                mode="markers+text", text=text, textposition="top center",
                name="PCA outlier candidate", hovertext=hover_text,
                hovertemplate="%{hovertext}<extra></extra>",
                marker={
                    "size": 9, "color": "#E15759", "symbol": "circle-open",
                    "line": {"color": "#E15759", "width": 2},
                },
            )
        )
    figure.write_html(
        output_path,
        include_plotlyjs=include_plotlyjs,
        full_html=True,
    )


def write_pca_not_run_outputs(
    distribution_root: Path,
    *,
    title: str,
    reason: str,
) -> None:
    """Write explicit placeholders when a valid PCA cannot be estimated."""
    message = f"PCA not run: {reason}"
    for filename, suffix in (
        ("pca_metadata.png", "PCA metadata"),
        ("pca_3d.png", "PCA 3D"),
    ):
        figure, axis = plt.subplots(figsize=(9, 5), facecolor="white")
        axis.axis("off")
        axis.text(0.5, 0.55, message, ha="center", va="center", fontsize=14)
        axis.set_title(f"{title} {suffix}")
        figure.tight_layout()
        figure.savefig(distribution_root / filename, dpi=220, bbox_inches="tight", facecolor="white")
        plt.close(figure)
    safe_title = html.escape(title)
    safe_message = html.escape(message)
    (distribution_root / "pca_3d_interactive.html").write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>{safe_title} PCA</title></head><body>"
        f"<h2>{safe_title} PCA</h2><p>{safe_message}</p></body></html>",
        encoding="utf-8",
    )


def write_qc_outputs(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    output_root: Path,
    *,
    title: str,
    include_descriptive: bool = True,
    include_plotlyjs: bool | str = True,
) -> QCOutputs:
    descriptive_root = output_root / "01_descriptive_statistics"
    distribution_root = (
        output_root / "02_data_distribution" if include_descriptive else output_root
    )
    if include_descriptive:
        descriptive_root.mkdir(parents=True, exist_ok=True)
    distribution_root.mkdir(parents=True, exist_ok=True)

    sample_statistics = calculate_sample_statistics(matrix, metadata)
    pca_status = "completed"
    try:
        pca_scores, pca_features, pca_removed = calculate_pca(matrix, metadata)
    except ValueError as error:
        pca_status = f"not_run: {error}"
        values = matrix.to_numpy(dtype=float)
        complete = np.isfinite(values).all(axis=0)
        complete_values = values[:, complete]
        nonconstant = (
            np.nanvar(complete_values, axis=0, ddof=0) > 0
            if complete_values.shape[1]
            else np.zeros(0, dtype=bool)
        )
        pca_features = int(nonconstant.sum())
        pca_removed = int((~complete).sum() + (~nonconstant).sum())
        pca_scores = metadata.copy()
        for component in ("PC1", "PC2", "PC3"):
            pca_scores[component] = np.nan
            pca_scores[f"{component}_explained_variance_ratio"] = np.nan
        pca_scores["pca_distance_pc1_pc3"] = np.nan
        pca_scores["pca_distance_robust_z"] = np.nan
        pca_scores["pca_outlier_candidate"] = False
        pca_scores["pca_status"] = pca_status
    sample_statistics["pca_distance_pc1_pc3"] = pca_scores["pca_distance_pc1_pc3"]
    sample_statistics["pca_distance_robust_z"] = pca_scores["pca_distance_robust_z"]
    sample_statistics["pca_outlier_candidate"] = pca_scores["pca_outlier_candidate"]
    sample_statistics["outlier_candidate"] = (
        sample_statistics["descriptive_outlier_candidate"]
        | sample_statistics["pca_outlier_candidate"]
    )
    pca_scores["outlier_candidate"] = sample_statistics["outlier_candidate"]

    paths: list[Path] = []
    gene_statistics: pd.DataFrame | None = None
    if include_descriptive:
        gene_statistics = calculate_gene_statistics(matrix)
        summary = pd.DataFrame([
            {"metric": "samples", "value": matrix.shape[0]},
            {"metric": "genes", "value": matrix.shape[1]},
            {"metric": "missing_cells", "value": int(matrix.isna().sum().sum())},
            {"metric": "missing_fraction", "value": float(matrix.isna().mean().mean())},
            {"metric": "pca_complete_nonconstant_genes", "value": pca_features},
            {"metric": "pca_genes_excluded_missing_or_constant", "value": pca_removed},
            {"metric": "pca_status", "value": pca_status},
            {"metric": "outlier_candidates", "value": int(sample_statistics["outlier_candidate"].sum())},
            {"metric": "expression_scale", "value": "stored values; no additional transform"},
        ])
        summary.to_csv(descriptive_root / "data_summary.csv", index=False)
        sample_statistics.reset_index(drop=True).to_csv(
            descriptive_root / "sample_statistics.csv", index=False
        )
        gene_statistics.to_csv(descriptive_root / "gene_statistics.csv", index=False)
        sample_statistics.loc[sample_statistics["outlier_candidate"]].reset_index(drop=True).to_csv(
            descriptive_root / "outlier_candidates.csv", index=False
        )
        plot_intensity_distribution(
            matrix, sample_statistics,
            descriptive_root / "intensity_distribution.png", title,
        )
        paths.extend([
            descriptive_root / "data_summary.csv",
            descriptive_root / "intensity_distribution.png",
            descriptive_root / "sample_statistics.csv",
            descriptive_root / "gene_statistics.csv",
            descriptive_root / "outlier_candidates.csv",
        ])

    pca_scores.reset_index(drop=True).to_csv(distribution_root / "pca_scores.csv", index=False)
    plot_heatmap(matrix, sample_statistics, distribution_root / "heatmap.png", title)
    if pca_status == "completed":
        plot_pca_metadata(pca_scores, distribution_root / "pca_metadata.png", f"{title} PCA metadata")
        plot_pca_3d(pca_scores, distribution_root / "pca_3d.png", f"{title} PCA")
        plot_pca_3d_interactive(
            pca_scores,
            distribution_root / "pca_3d_interactive.html",
            f"{title} PCA (interactive)",
            include_plotlyjs=include_plotlyjs,
        )
    else:
        write_pca_not_run_outputs(
            distribution_root,
            title=title,
            reason=pca_status.removeprefix("not_run: "),
        )
    paths.extend([
        distribution_root / "heatmap.png",
        distribution_root / "pca_metadata.png",
        distribution_root / "pca_3d.png",
        distribution_root / "pca_3d_interactive.html",
        distribution_root / "pca_scores.csv",
    ])
    return QCOutputs(
        sample_statistics=sample_statistics,
        gene_statistics=gene_statistics,
        pca_scores=pca_scores,
        pca_features_used=pca_features,
        paths=tuple(paths),
    )

