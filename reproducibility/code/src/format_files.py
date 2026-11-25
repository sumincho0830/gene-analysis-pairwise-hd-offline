"""Build sample-by-feature matrices for the first preprocessing stage.

This module contains only file-format conversion. It does not filter samples,
filter genes, transform values, or run PyWGCNA.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


RUN_ROOT = Path(__file__).resolve().parents[2]

PROTEOMICS_RAW_DIR = RUN_ROOT / "data/proteomics/cohorts/1_raw"
PROTEOMICS_CORRECTED_DIR = (
    RUN_ROOT / "data/proteomics/cohorts/2_batch_corrected"
)
TRANSCRIPTOMICS_DIR = RUN_ROOT / "data/transcriptomics/cohorts"
CELL_TYPES_FILE = RUN_ROOT / "data/metadata/cell_types.csv"

FORMATTED_OUTPUT_DIR = RUN_ROOT / "formatted"
PROTEOMICS_RAW_OUTPUT_DIR = FORMATTED_OUTPUT_DIR / "00_raw/proteomics/pywgcna_input"
PROTEOMICS_CORRECTED_OUTPUT_DIR = (
    FORMATTED_OUTPUT_DIR / "01_batch_corrected/proteomics/pywgcna_input"
)
TRANSCRIPTOMICS_OUTPUT_DIR = (
    FORMATTED_OUTPUT_DIR / "00_raw/transcriptomics/pywgcna_input"
)
TRANSCRIPTOMICS_CORRECTED_DIR = RUN_ROOT / "data/Transcriptomics_Corrected"
TRANSCRIPTOMICS_CORRECTED_OUTPUT_DIR = (
    FORMATTED_OUTPUT_DIR / "01_batch_corrected/transcriptomics"
)

PROTEOMICS_ANNOTATION_COLUMNS = (
    "ProteinColumn",
    "PROBEID",
    "SYMBOL",
    "UNIPROT",
)
SPECIMEN_PATTERN = re.compile(r"_(?P<specimen>CSF|PB)(?:\d+)?(?:_|$)")
TRANSCRIPTOMICS_COLUMN_ALIASES = {
    "Classical Monocytes": "CD14 Mono",
    "Non-classical Monocytes": "CD16 Mono",
}


def clean_name(name: object) -> str:
    """Return the historical file-safe cell-type name."""

    return str(name).replace(" ", "_").replace("-", "_").replace("/", "_")


def extract_specimen(sample_id: object) -> str:
    """Extract CSF or PB from the end of a sample identifier."""

    match = SPECIMEN_PATTERN.search(str(sample_id))
    if match is None:
        raise ValueError(f"Cannot determine specimen from sample ID: {sample_id!r}")
    return match.group("specimen")


def _require_columns(
    frame: pd.DataFrame,
    required: Iterable[str],
    *,
    source: str,
) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{source}: missing required columns {missing}")


def _validate_feature_ids(feature_ids: pd.Index, *, source: str) -> None:
    text = feature_ids.astype(str)
    if (text.str.strip() == "").any() or (text.str.lower() == "nan").any():
        raise ValueError(f"{source}: feature IDs contain an empty value")
    if text.duplicated().any():
        duplicated = text[text.duplicated()].unique().tolist()[:5]
        raise ValueError(f"{source}: duplicated feature IDs {duplicated}")


def _as_numeric(
    frame: pd.DataFrame,
    *,
    source: str,
    allow_missing: bool,
) -> pd.DataFrame:
    try:
        numeric = frame.apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{source}: expression values must be numeric") from exc

    values = numeric.to_numpy(dtype=float, copy=False)
    if np.isinf(values).any():
        raise ValueError(f"{source}: expression values contain infinity")
    if not allow_missing and np.isnan(values).any():
        raise ValueError(f"{source}: expression values contain NaN")
    return numeric


def format_proteomics_dataframe(
    frame: pd.DataFrame,
    *,
    source: str = "<proteomics dataframe>",
    allow_missing: bool = False,
) -> pd.DataFrame:
    """Convert one cohort table to samples x PROBEIDs.

    The first source row stores original sample identifiers and is excluded.
    No log transform, normalization, filtering, or aggregation is performed.
    ``allow_missing=True`` preserves raw-platform NaN values without imputation.
    """

    _require_columns(frame, PROTEOMICS_ANNOTATION_COLUMNS, source=source)
    if len(frame) < 2:
        raise ValueError(f"{source}: expected one metadata row and feature rows")

    feature_rows = frame.iloc[1:]
    feature_ids = pd.Index(feature_rows["PROBEID"].astype(str), name=None)
    _validate_feature_ids(feature_ids, source=source)

    sample_columns = [
        column
        for column in frame.columns
        if column not in PROTEOMICS_ANNOTATION_COLUMNS
    ]
    if not sample_columns:
        raise ValueError(f"{source}: no sample columns found")

    expression = feature_rows.loc[:, sample_columns].T
    expression.columns = feature_ids
    expression.columns.name = None
    expression.index = expression.index.astype(str)
    expression.index.name = "PROBEID"  # Preserve the existing output contract.
    expression = _as_numeric(
        expression,
        source=source,
        allow_missing=allow_missing,
    )
    return expression


def format_transcriptomics_sample(
    frame: pd.DataFrame,
    *,
    sample_id: str,
    cell_type: str,
    source: str = "<transcriptomics dataframe>",
) -> pd.DataFrame:
    """Convert one sample/cell-type column to a one-row samples x genes table."""

    source_column = TRANSCRIPTOMICS_COLUMN_ALIASES.get(cell_type, cell_type)
    _require_columns(frame, ("gene_symbol", source_column), source=source)
    gene_ids = pd.Index(frame["gene_symbol"].astype(str), name=None)
    _validate_feature_ids(gene_ids, source=source)

    expression = pd.DataFrame(
        [frame[source_column].to_numpy()],
        index=pd.Index([str(sample_id)], name="gene_symbol"),
        columns=gene_ids,
    )
    expression.columns.name = None
    return _as_numeric(expression, source=source, allow_missing=False)


def split_by_specimen(expression: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Split a samples x features table without silently dropping sample IDs."""

    specimens = pd.Series(
        [extract_specimen(sample_id) for sample_id in expression.index],
        index=expression.index,
    )
    return {
        "ALL": expression.copy(),
        "CSF": expression.loc[specimens.eq("CSF")].copy(),
        "PB": expression.loc[specimens.eq("PB")].copy(),
    }


def _concat_expression_frames(
    frames: list[pd.DataFrame],
    *,
    source: str,
    require_same_features: bool,
) -> pd.DataFrame:
    if not frames:
        raise ValueError(f"{source}: no expression tables were produced")

    if require_same_features:
        expected = frames[0].columns
        for index, frame in enumerate(frames[1:], start=2):
            if not frame.columns.equals(expected):
                raise ValueError(
                    f"{source}: table {index} has different feature IDs or order"
                )

    combined = pd.concat(frames, axis=0, sort=False)
    if require_same_features and combined.isna().any().any():
        raise ValueError(
            f"{source}: missing values appeared while combining expression tables"
        )
    if combined.index.duplicated().any():
        duplicated = combined.index[combined.index.duplicated()].unique().tolist()[:5]
        raise ValueError(f"{source}: duplicated sample IDs {duplicated}")
    return combined


def load_proteomics_matrix(
    root: Path,
    *,
    allow_missing: bool = False,
) -> pd.DataFrame:
    """Read all proteomics cohort CSVs once and combine their sample rows."""

    files = sorted(Path(root).glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {root}")
    frames = [
        format_proteomics_dataframe(
            pd.read_csv(path),
            source=str(path),
            allow_missing=allow_missing,
        )
        for path in files
    ]
    # Raw platforms do not all measure the same probes. pandas aligns them by
    # PROBEID and leaves unmeasured platform-specific probes as NaN, matching
    # the existing formatted-file contract.
    return _concat_expression_frames(
        frames,
        source=str(root),
        require_same_features=False,
    )


def read_cell_types(path: Path = CELL_TYPES_FILE) -> list[str]:
    """Read the ordered cell-type list used by transcriptomics inputs."""

    values = pd.read_csv(path, header=None).iloc[:, 0].dropna().astype(str).tolist()
    if not values:
        raise ValueError(f"{path}: no cell types found")
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicated cell types")
    return values


def load_transcriptomics_matrices(
    root: Path,
    *,
    cell_types: Iterable[str],
) -> dict[str, pd.DataFrame]:
    """Read each transcriptomics sample once and build one matrix per cell type."""

    files = sorted(Path(root).glob("*/*.csv"))
    if not files:
        raise FileNotFoundError(f"No cohort CSV files found in {root}")

    ordered_cell_types = list(cell_types)
    rows: dict[str, list[pd.DataFrame]] = {
        cell_type: [] for cell_type in ordered_cell_types
    }
    for path in files:
        frame = pd.read_csv(path)
        for cell_type in ordered_cell_types:
            source_column = TRANSCRIPTOMICS_COLUMN_ALIASES.get(cell_type, cell_type)
            if source_column not in frame.columns:
                continue
            rows[cell_type].append(
                format_transcriptomics_sample(
                    frame,
                    sample_id=path.stem,
                    cell_type=cell_type,
                    source=str(path),
                )
            )

    matrices: dict[str, pd.DataFrame] = {}
    for cell_type, cell_rows in rows.items():
        matrices[cell_type] = _concat_expression_frames(
            cell_rows,
            source=f"{root} [{cell_type}]",
            require_same_features=True,
        )
    return matrices


def _write_matrix(expression: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    expression.to_csv(path)
    return path


def build_proteomics_expression_matrices(
    *,
    root: Path,
    output_dir: Path,
    allow_missing: bool = False,
) -> dict[str, Path]:
    """Write ALL, CSF, and PB proteomics matrices without value imputation."""

    split = split_by_specimen(
        load_proteomics_matrix(root, allow_missing=allow_missing)
    )
    names = {
        "ALL": "PROT_ALL_expression_matrix.csv",
        "CSF": "PROT_CSF_expression_matrix.csv",
        "PB": "PROT_PB_expression_matrix.csv",
    }
    return {
        specimen: _write_matrix(frame, Path(output_dir) / specimen / names[specimen])
        for specimen, frame in split.items()
    }


def build_transcriptomics_expression_matrices(
    *,
    root: Path = TRANSCRIPTOMICS_DIR,
    output_dir: Path = TRANSCRIPTOMICS_OUTPUT_DIR,
    cell_types_file: Path = CELL_TYPES_FILE,
) -> dict[str, dict[str, Path]]:
    """Write ALL, CSF, and PB matrices for every configured cell type."""

    matrices = load_transcriptomics_matrices(
        root,
        cell_types=read_cell_types(cell_types_file),
    )
    written: dict[str, dict[str, Path]] = {}
    for cell_type, expression in matrices.items():
        safe_name = clean_name(cell_type)
        split = split_by_specimen(expression)
        names = {
            "ALL": f"{safe_name}_expression_matrix.csv",
            "CSF": f"{safe_name}_csf_expression_matrix.csv",
            "PB": f"{safe_name}_pb_expression_matrix.csv",
        }
        written[cell_type] = {
            specimen: _write_matrix(
                frame,
                Path(output_dir) / specimen / names[specimen],
            )
            for specimen, frame in split.items()
        }
    return written


def build_corrected_transcriptomics_expression_matrices(
    *,
    corrected_root: Path = TRANSCRIPTOMICS_CORRECTED_DIR,
    output_dir: Path = TRANSCRIPTOMICS_CORRECTED_OUTPUT_DIR,
    cell_types_file: Path = CELL_TYPES_FILE,
) -> dict[str, dict[str, dict[str, Path]]]:
    """Format every corrected transcriptomics threshold without mixing them.

    Each threshold is read only from its ``Batch_Corrected_CSVs`` directory and
    written under its own directory. Therefore the threshold used to retain a
    cell type remains explicit in every downstream matrix path.
    """

    threshold_dirs = sorted(
        path
        for path in Path(corrected_root).glob("Threshold_*_Cells")
        if path.is_dir()
    )
    if not threshold_dirs:
        raise FileNotFoundError(
            f"No Threshold_*_Cells directories found in {corrected_root}"
        )

    written: dict[str, dict[str, dict[str, Path]]] = {}
    for threshold_dir in threshold_dirs:
        source_dir = threshold_dir / "Batch_Corrected_CSVs"
        if not source_dir.is_dir():
            raise FileNotFoundError(f"Missing batch-corrected input: {source_dir}")
        output_name = clean_name(threshold_dir.name).lower()
        written[threshold_dir.name] = build_transcriptomics_expression_matrices(
            root=source_dir,
            output_dir=Path(output_dir) / output_name / "pywgcna_input",
            cell_types_file=cell_types_file,
        )
    return written


def format_all_files() -> dict[str, object]:
    """Run only stage 1 and return the paths of all written matrices."""

    return {
        "proteomics_raw": build_proteomics_expression_matrices(
            root=PROTEOMICS_RAW_DIR,
            output_dir=PROTEOMICS_RAW_OUTPUT_DIR,
            allow_missing=True,
        ),
        "proteomics_batch_corrected": build_proteomics_expression_matrices(
            root=PROTEOMICS_CORRECTED_DIR,
            output_dir=PROTEOMICS_CORRECTED_OUTPUT_DIR,
        ),
        "transcriptomics": build_transcriptomics_expression_matrices(),
        "transcriptomics_batch_corrected": (
            build_corrected_transcriptomics_expression_matrices()
        ),
    }


if __name__ == "__main__":
    written = format_all_files()
    print(f"Formatted data written under {FORMATTED_OUTPUT_DIR}")
