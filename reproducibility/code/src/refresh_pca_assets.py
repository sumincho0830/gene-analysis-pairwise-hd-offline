"""Refresh PCA image assets in an existing offline dashboard release.

The release uses content-addressed lossless WebP files.  This script matches
records and artifact IDs against the source dashboard, converts the current
PNG figures, updates only the affected artifact paths in ``dashboard-data.js``,
and rebuilds the release checksum list.  Other dashboard UI/data are preserved.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any


DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
DEFAULT_ARTIFACT_IDS = ("before_pca", "after_pca")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_bundle(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8").strip()
    if not text.startswith(DATA_PREFIX) or not text.endswith(";"):
        raise ValueError(f"Unexpected dashboard data wrapper: {path}")
    return json.loads(text[len(DATA_PREFIX) : -1])


def write_bundle(path: Path, bundle: dict[str, Any]) -> None:
    path.write_text(
        DATA_PREFIX
        + json.dumps(bundle, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )


def convert_lossless_webp(task: tuple[str, str]) -> dict[str, Any]:
    from PIL import Image

    source = Path(task[0])
    target = Path(task[1])
    Image.MAX_IMAGE_PIXELS = None
    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as image:
        mode = "RGBA" if "A" in image.getbands() else "RGB"
        pixels = image.convert(mode)
        expected_hash = hashlib.sha256(pixels.tobytes()).hexdigest()
        expected_size = pixels.size
        pixels.save(target, "WEBP", lossless=True, method=2, exact=True)
    with Image.open(target) as converted:
        decoded = converted.convert(mode)
        actual_hash = hashlib.sha256(decoded.tobytes()).hexdigest()
        if decoded.size != expected_size or actual_hash != expected_hash:
            raise ValueError(f"Lossless pixel validation failed: {source} -> {target}")
    return {"source": str(source), "target": str(target)}


def rebuild_checksums(release_root: Path) -> int:
    paths = sorted(
        path
        for path in release_root.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS.txt"
        and ".git" not in path.relative_to(release_root).parts
        and "__pycache__" not in path.relative_to(release_root).parts
        and path.suffix.lower() != ".pyc"
    )
    (release_root / "SHA256SUMS.txt").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(release_root).as_posix()}"
            for path in paths
        )
        + "\n",
        encoding="utf-8",
    )
    return len(paths) + 1


def refresh(
    source_run: Path,
    release_root: Path,
    artifact_ids: set[str],
    workers: int,
) -> dict[str, Any]:
    source_dashboard = source_run / "dashboard"
    source_bundle = read_bundle(source_dashboard / "dashboard-data.js")
    release_data_path = release_root / "dashboard-data.js"
    release_bundle = read_bundle(release_data_path)
    source_records = {record["id"]: record for record in source_bundle["records"]}

    conversions: dict[str, tuple[Path, Path]] = {}
    updated_references = 0
    missing: list[str] = []
    for release_record in release_bundle["records"]:
        source_record = source_records.get(release_record["id"])
        if source_record is None:
            missing.append(f"record:{release_record['id']}")
            continue
        source_artifacts = {
            artifact["id"]: artifact for artifact in source_record.get("artifacts", [])
        }
        for release_artifact in release_record.get("artifacts", []):
            artifact_id = release_artifact.get("id")
            if artifact_id not in artifact_ids or release_artifact.get("kind") != "image":
                continue
            source_artifact = source_artifacts.get(artifact_id)
            if source_artifact is None:
                missing.append(f"artifact:{release_record['id']}:{artifact_id}")
                continue
            source = (source_dashboard / source_artifact["path"]).resolve()
            if not source.is_file():
                missing.append(f"file:{source}")
                continue
            digest = sha256(source)[:16]
            target = release_root / "assets" / artifact_id / f"{digest}.webp"
            target.resolve().relative_to(release_root.resolve())
            conversions[target.as_posix().lower()] = (source, target)
            new_path = target.relative_to(release_root).as_posix()
            if release_artifact.get("path") != new_path:
                release_artifact["path"] = new_path
                updated_references += 1

    tasks = [
        (str(source), str(target))
        for source, target in conversions.values()
        if not target.is_file()
    ]
    converted = 0
    with ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = [executor.submit(convert_lossless_webp, task) for task in tasks]
        for future in as_completed(futures):
            future.result()
            converted += 1
            if converted % 25 == 0 or converted == len(futures):
                print(f"Converted {converted}/{len(futures)} PCA assets", flush=True)

    now = datetime.now(timezone.utc).isoformat()
    release_bundle["generated_at_utc"] = now
    write_bundle(release_data_path, release_bundle)

    manifest_path = release_root / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["last_pca_asset_refresh_utc"] = now
    manifest["pca_asset_references_refreshed"] = updated_references
    manifest["pca_asset_files_converted"] = converted
    manifest["pca_asset_missing_sources"] = len(missing)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    manifest["release_files"] = rebuild_checksums(release_root)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    rebuild_checksums(release_root)

    return {
        "artifact_ids": sorted(artifact_ids),
        "updated_references": updated_references,
        "unique_assets": len(conversions),
        "converted_assets": converted,
        "missing_sources": missing,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument(
        "--release-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
    )
    parser.add_argument(
        "--artifact-id",
        action="append",
        dest="artifact_ids",
        help="Artifact ID to refresh; repeat for multiple IDs.",
    )
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = refresh(
        args.source_run.resolve(),
        args.release_root.resolve(),
        set(args.artifact_ids or DEFAULT_ARTIFACT_IDS),
        args.workers,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
