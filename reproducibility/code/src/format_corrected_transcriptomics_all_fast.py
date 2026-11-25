"""Format corrected transcriptomics inputs into ALL sample-by-gene matrices.

This specialized formatter reads every cohort CSV once per threshold. It keeps
the existing gene order and expression values, maps source cell-type aliases,
and normalizes source specimen names (BLD -> PB) in sample identifiers.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from format_files import (
    CELL_TYPES_FILE,
    TRANSCRIPTOMICS_CORRECTED_OUTPUT_DIR,
    clean_name,
    extract_specimen,
    read_cell_types,
)


SOURCE_ALIASES = {
    "Classical Monocytes": "CD14 Mono",
    "Non-classical Monocytes": "CD16 Mono",
}
SAMPLE_PATTERN = re.compile(
    r"^bulk_counts_(?P<base>.+)_(?P<specimen>CSF\d*|PB|BLD)_(?P<diagnosis>[^_]+)$"
)


def normalized_sample_id(path: Path) -> tuple[str, str]:
    match = SAMPLE_PATTERN.match(path.stem)
    if match is None:
        raise ValueError(f"Unrecognized corrected transcriptomics filename: {path.name}")
    source_specimen = match.group("specimen")
    specimen = "CSF" if source_specimen.startswith("CSF") else "PB"
    output_suffix = source_specimen if specimen == "CSF" else "PB"
    return f"{match.group('base')}_{output_suffix}", specimen


def format_threshold(
    threshold_dir: Path,
    cell_types: list[str],
    *,
    output_root: Path,
) -> list[Path]:
    source_dir = threshold_dir / "Batch_Corrected_CSVs"
    paths = sorted(source_dir.glob("*/*.csv"))
    if not paths:
        raise FileNotFoundError(f"No batch-corrected CSVs in {source_dir}")

    sample_ids: dict[str, list[str]] = {cell_type: [] for cell_type in cell_types}
    rows: dict[str, list[np.ndarray]] = {cell_type: [] for cell_type in cell_types}
    expected_genes: pd.Index | None = None

    for path in paths:
        frame = pd.read_csv(path)
        if "gene_symbol" not in frame.columns:
            raise ValueError(f"{path}: missing gene_symbol")
        genes = pd.Index(frame["gene_symbol"].astype(str))
        if expected_genes is None:
            expected_genes = genes
        elif not genes.equals(expected_genes):
            raise ValueError(f"{path}: gene order differs from the first input")

        sample_id, _ = normalized_sample_id(path)
        for cell_type in cell_types:
            source_column = SOURCE_ALIASES.get(cell_type, cell_type)
            if source_column not in frame.columns:
                continue
            values = pd.to_numeric(frame[source_column], errors="raise").to_numpy(
                dtype=float, copy=True
            )
            if not np.isfinite(values).all():
                raise ValueError(f"{path}: non-finite values in {source_column}")
            sample_ids[cell_type].append(sample_id)
            rows[cell_type].append(values)

    assert expected_genes is not None
    threshold_name = threshold_dir.name.lower()
    base_output_dir = output_root / threshold_name / "pywgcna_input"
    written: list[Path] = []
    for cell_type in cell_types:
        ids = sample_ids[cell_type]
        if not ids:
            raise ValueError(f"{threshold_dir.name}: no rows for {cell_type}")
        if len(ids) != len(set(ids)):
            raise ValueError(f"{threshold_dir.name}: duplicate sample IDs for {cell_type}")
        matrix = pd.DataFrame(
            np.vstack(rows[cell_type]),
            index=pd.Index(ids, name="gene_symbol"),
            columns=expected_genes,
        )
        for specimen in ("ALL", "CSF", "PB"):
            subset = (
                matrix
                if specimen == "ALL"
                else matrix.loc[
                    [extract_specimen(sample_id) == specimen for sample_id in matrix.index]
                ]
            )
            if subset.empty:
                raise ValueError(
                    f"{threshold_dir.name}: no {specimen} rows for {cell_type}"
                )
            output_dir = base_output_dir / specimen
            output_dir.mkdir(parents=True, exist_ok=True)
            output = output_dir / f"{clean_name(cell_type)}_expression_matrix.csv"
            subset.to_csv(output)
            written.append(output)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "data/Transcriptomics_Corrected",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=TRANSCRIPTOMICS_CORRECTED_OUTPUT_DIR,
    )
    parser.add_argument("--cell-types-file", type=Path, default=CELL_TYPES_FILE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_root = args.input_root.resolve()
    output_root = args.output_root.resolve()
    cell_types = read_cell_types(args.cell_types_file.resolve())
    written: list[Path] = []
    for threshold_dir in sorted(input_root.glob("Threshold_*_Cells")):
        written.extend(
            format_threshold(threshold_dir, cell_types, output_root=output_root)
        )
        print(f"Finished {threshold_dir.name}", flush=True)
    print(f"Saved {len(written)} ALL/CSF/PB matrices under {output_root}")


if __name__ == "__main__":
    main()
