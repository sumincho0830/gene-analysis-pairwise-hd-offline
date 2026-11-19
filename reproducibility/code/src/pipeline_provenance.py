"""Small, reusable provenance writer for dated pipeline runs."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
from typing import Any, Iterable


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _file_record(path: Path, run_root: Path) -> dict[str, Any]:
    resolved = path.resolve()
    try:
        display_path = resolved.relative_to(run_root.resolve()).as_posix()
    except ValueError:
        display_path = str(resolved)
    return {
        "path": display_path,
        "size_bytes": resolved.stat().st_size,
        "sha256": _sha256(resolved),
    }


def _records(paths: Iterable[Path], run_root: Path) -> list[dict[str, Any]]:
    return [_file_record(path, run_root) for path in paths if path.exists()]


def _package_versions(packages: Iterable[str]) -> dict[str, str]:
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def write_run_manifest(
    output_path: Path,
    *,
    pipeline_name: str,
    run_root: Path,
    parameters: dict[str, Any],
    generation_rules: dict[str, Any],
    input_paths: Iterable[Path],
    output_paths: Iterable[Path],
    code_paths: Iterable[Path],
    packages: Iterable[str] = ("numpy", "pandas", "matplotlib", "plotly"),
) -> None:
    """Record enough state to audit whether a later run used the same inputs and rules."""
    payload = {
        "schema_version": 1,
        "pipeline": pipeline_name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(run_root.resolve()),
        "parameters": parameters,
        "generation_rules": generation_rules,
        "inputs": _records(input_paths, run_root),
        "outputs": _records(output_paths, run_root),
        "code": _records(code_paths, run_root),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": _package_versions(packages),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
