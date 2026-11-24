"""Build the publication-style static version of the final CSF Jaccard heatmap.

The figure uses the same saved overlap index, condition, cell-type order,
protein-module order, and within-cell RNA-module ordering as final-jaccard.html.
No similarity values are recalculated.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
INDEX_PATH = ROOT / "cross-omics-overlap" / "overlap-index.js"
OUTPUT_DIR = ROOT / "figures" / "final_jaccard"
PREFIX = "window.CROSS_OMICS_OVERLAP_INDEX="

THRESHOLD = "threshold_5_cells"
STRATEGY = "02_hvg_top9000"
REMOVAL = "after"
SPECIMEN = "CSF"
MODE = "cell_separate"

CELL_ORDER = [
    "B_cell",
    "CD4_T_cell",
    "CD8_T_cell",
    "cDCs",
    "Classical_Monocytes",
    "Macrophages",
    "MAIT_gDT_cell",
    "NK_cell",
    "Non_classical_Monocytes",
    "pDCs",
    "Plasmablast",
    "Proliferating_T_cell",
]

FILE_STEM = (
    "CSF_threshold_5_cells_hvg_top9000_after_outlier_"
    "kME_weighted_Jaccard"
)


def read_index(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8").strip()
    if not text.startswith(PREFIX) or not text.endswith(";"):
        raise ValueError(f"Unexpected overlap-index wrapper: {path}")
    return json.loads(text[len(PREFIX) : -1])


def dimensions(record_id: str, entry: dict[str, Any]) -> dict[str, str]:
    parts = record_id.split("::")
    after = parts[-1] == "after_outlier"
    base = parts[:-1] if after else parts
    return {
        "threshold": str(entry.get("threshold") or base[1]),
        "strategy": str(entry.get("strategy") or base[2]),
        "mode": str(entry.get("mode") or base[3]),
        "specimen": str(entry.get("specimen") or base[4]),
        "cell_type": str(entry.get("cell_type") or base[5]),
        "removal": "after" if after else "before",
    }


def pair_objects(entry: dict[str, Any]) -> list[dict[str, Any]]:
    columns = entry["pair_columns"]
    return [dict(zip(columns, row)) for row in entry["pairs"]]


def selected_records(index: dict[str, Any]) -> list[tuple[str, dict[str, Any], dict[str, str]]]:
    order = {cell: position for position, cell in enumerate(CELL_ORDER)}
    selected: list[tuple[str, dict[str, Any], dict[str, str]]] = []
    for record_id, entry in index["records"].items():
        dims = dimensions(record_id, entry)
        if (
            dims["specimen"] == SPECIMEN
            and dims["mode"] == MODE
            and dims["cell_type"].lower() != "other"
            and dims["threshold"] == THRESHOLD
            and dims["strategy"] == STRATEGY
            and dims["removal"] == REMOVAL
        ):
            selected.append((record_id, entry, dims))
    selected.sort(key=lambda item: (order.get(item[2]["cell_type"], 999), item[2]["cell_type"]))
    missing = [cell for cell in CELL_ORDER if cell not in {item[2]["cell_type"] for item in selected}]
    if missing:
        raise RuntimeError(f"Missing final-condition cell types: {missing}")
    return selected


def prepare_matrix(index: dict[str, Any]) -> dict[str, Any]:
    records = selected_records(index)
    protein_modules = sorted(
        {module for _, entry, _ in records for module in entry["proteomics_modules"]}
    )
    rows: list[dict[str, str]] = []
    pair_lookup: dict[tuple[str, str, str], dict[str, Any]] = {}

    for record_id, entry, dims in records:
        pairs = pair_objects(entry)
        for pair in pairs:
            pair_lookup[(record_id, pair["transcript_module"], pair["proteomics_module"])] = pair
        module_max = {
            module: max(
                float(pair_lookup.get((record_id, module, protein), {}).get("weighted_jaccard", 0) or 0)
                for protein in protein_modules
            )
            for module in entry["transcript_modules"]
        }
        modules = sorted(entry["transcript_modules"], key=lambda module: (-module_max[module], module))
        for module in modules:
            rows.append(
                {
                    "record_id": record_id,
                    "cell_type": dims["cell_type"],
                    "transcript_module": module,
                }
            )

    matrix = np.full((len(rows), len(protein_modules)), np.nan, dtype=float)
    transcript_best = np.zeros(matrix.shape, dtype=bool)
    reciprocal_best = np.zeros(matrix.shape, dtype=bool)
    for y, row in enumerate(rows):
        for x, protein_module in enumerate(protein_modules):
            pair = pair_lookup.get(
                (row["record_id"], row["transcript_module"], protein_module)
            )
            if pair is None:
                continue
            matrix[y, x] = float(pair["weighted_jaccard"])
            transcript_best[y, x] = bool(pair.get("transcript_best"))
            reciprocal_best[y, x] = bool(pair.get("reciprocal_best"))

    if not np.isfinite(matrix).any():
        raise RuntimeError("The final Jaccard matrix has no finite values")
    return {
        "records": records,
        "rows": rows,
        "protein_modules": protein_modules,
        "matrix": matrix,
        "transcript_best": transcript_best,
        "reciprocal_best": reciprocal_best,
    }


def pretty(value: str) -> str:
    return value.replace("_", " ")


def plot(data: dict[str, Any], output_dir: Path) -> list[Path]:
    rows = data["rows"]
    protein_modules = data["protein_modules"]
    matrix = data["matrix"]
    transcript_best = data["transcript_best"]
    reciprocal_best = data["reciprocal_best"]
    n_rows, n_columns = matrix.shape
    vmax = float(np.nanmax(matrix))

    cmap = LinearSegmentedColormap.from_list(
        "final_jaccard", ["#ffffff", "#d7e7ee", "#8fb8c9", "#1f6485", "#153f56"]
    )
    cmap.set_bad("#eeeeee")
    cell_ranges: list[tuple[str, int, int]] = []
    start = 0
    for position in range(1, n_rows + 1):
        if position == n_rows or rows[position]["cell_type"] != rows[start]["cell_type"]:
            cell_ranges.append((rows[start]["cell_type"], start, position))
            start = position
    possible_splits = [last for _, _, last in cell_ranges[:-1]]
    split = min(possible_splits, key=lambda value: abs(value - n_rows / 2))
    panels = [(0, split), (split, n_rows)]
    figure, axes = plt.subplots(1, 2, figsize=(18.5, 20.5), facecolor="white")
    mesh = None
    for panel_index, (axis, (first_row, last_row)) in enumerate(zip(axes, panels)):
        panel_matrix = matrix[first_row:last_row]
        mesh = axis.pcolormesh(
            np.ma.masked_invalid(panel_matrix),
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            shading="flat",
            edgecolors="#edf0f2",
            linewidth=0.16,
            rasterized=False,
        )
        panel_rows = last_row - first_row
        axis.set_xlim(0, n_columns)
        axis.set_ylim(panel_rows, 0)
        axis.set_aspect("auto")
        axis.set_xticks(np.arange(n_columns) + 0.5)
        axis.set_xticklabels(protein_modules, rotation=55, ha="right", fontsize=7)
        axis.xaxis.tick_top()
        axis.xaxis.set_label_position("top")
        axis.set_xlabel("Proteomics module", fontsize=9, labelpad=10)
        axis.set_yticks(np.arange(panel_rows) + 0.5)
        axis.set_yticklabels(
            [row["transcript_module"] for row in rows[first_row:last_row]],
            fontsize=4.5,
        )
        axis.tick_params(axis="y", length=0, pad=2)
        axis.tick_params(axis="x", length=0, pad=3)
        axis.text(
            -0.19,
            1.015,
            chr(ord("A") + panel_index),
            transform=axis.transAxes,
            fontsize=12,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

        transform = blended_transform_factory(axis.transAxes, axis.transData)
        for cell_type, cell_first, cell_last in cell_ranges:
            if cell_first < first_row or cell_last > last_row:
                continue
            local_first = cell_first - first_row
            local_last = cell_last - first_row
            if local_first > 0:
                axis.axhline(local_first, color="#52636e", linewidth=0.85, zorder=4)
            axis.text(
                -0.155,
                (local_first + local_last) / 2,
                pretty(cell_type),
                transform=transform,
                ha="right",
                va="center",
                fontsize=7.0,
                fontweight="bold",
                color="#26343d",
                clip_on=False,
            )

        for global_y, x in np.argwhere(transcript_best[first_row:last_row]):
            y = int(global_y)
            source_y = first_row + y
            if reciprocal_best[source_y, x]:
                axis.add_patch(
                    Rectangle(
                        (x + 0.08, y + 0.08),
                        0.84,
                        0.84,
                        fill=False,
                        edgecolor="#f0a202",
                        linewidth=0.85,
                        zorder=6,
                    )
                )
            else:
                value = matrix[source_y, x]
                axis.plot(
                    x + 0.5,
                    y + 0.5,
                    marker="o",
                    markersize=1.8,
                    color="white" if value >= vmax * 0.55 else "#25323a",
                    markeredgewidth=0,
                    zorder=5,
                )
        for spine in axis.spines.values():
            spine.set_visible(False)

    if mesh is None:
        raise RuntimeError("No heatmap panel was drawn")
    colorbar_axis = figure.add_axes([0.38, 0.034, 0.24, 0.009])
    colorbar = figure.colorbar(mesh, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("|kME|-weighted Jaccard similarity", fontsize=8, labelpad=4)
    colorbar.ax.tick_params(labelsize=7, length=2)

    figure.text(
        0.07,
        0.992,
        "CSF transcriptome-proteome module similarity",
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
        color="#17232b",
    )
    figure.text(
        0.07,
        0.973,
        "Threshold: 5 cells  |  Features: HVG top 9,000  |  After removal of Ng_MOG1 and Ng_HD3\n"
        f"12 cell types (other excluded)  |  {n_rows} RNA modules x {n_columns} protein modules",
        ha="left",
        va="top",
        fontsize=9,
        color="#55636c",
        linespacing=1.35,
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#25323a",
            markeredgewidth=0,
            markersize=4,
            label="RNA-module best protein match",
        ),
        Patch(
            facecolor="none",
            edgecolor="#f0a202",
            linewidth=1.1,
            label="Reciprocal best match",
        ),
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.965, 0.986),
        frameon=False,
        fontsize=8,
        ncol=2,
        handlelength=1.3,
    )
    figure.text(
        0.5,
        0.006,
        "Descriptive similarity from saved module kME values; permutation p-values are not included.",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="#66717a",
    )
    figure.subplots_adjust(left=0.16, right=0.985, top=0.90, bottom=0.055, wspace=0.28)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir / f"{FILE_STEM}.{suffix}" for suffix in ("png", "pdf", "svg")]
    figure.savefig(paths[0], dpi=300, facecolor="white")
    figure.savefig(paths[1], facecolor="white")
    figure.savefig(paths[2], facecolor="white")
    plt.close(figure)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", type=Path, default=INDEX_PATH)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    data = prepare_matrix(read_index(args.index.resolve()))
    paths = plot(data, args.output_dir.resolve())
    matrix = data["matrix"]
    print(
        json.dumps(
            {
                "condition": {
                    "specimen": SPECIMEN,
                    "threshold": THRESHOLD,
                    "strategy": STRATEGY,
                    "outlier_removal": REMOVAL,
                },
                "cell_types": len(data["records"]),
                "rna_modules": int(matrix.shape[0]),
                "protein_modules": int(matrix.shape[1]),
                "module_pairs": int(np.isfinite(matrix).sum()),
                "maximum_weighted_jaccard": float(np.nanmax(matrix)),
                "outputs": [str(path) for path in paths],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
