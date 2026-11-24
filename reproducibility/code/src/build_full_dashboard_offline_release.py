"""Build a complete, cloneable offline bundle of the inspection dashboard.

Every artifact and downloadable file linked by ``dashboard-data.js`` is copied
into the release. Raster figures are converted to pixel-identical lossless
WebP to keep the Git repository practical. Plotly HTML files are rewritten to
use one bundled ``plotly.min.js`` instead of a CDN.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shutil
from typing import Any


DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
PLOTLY_CDN_RE = re.compile(
    r"<script\b[^>]*\bsrc=[\"']https://cdn\.plot\.ly/[^\"']+[\"'][^>]*></script>",
    flags=re.IGNORECASE,
)
CODE_FILES = (
    "code/src/build_experiment_inspection_dashboard.py",
    "code/src/build_full_dashboard_offline_release.py",
    "code/src/run_transcriptomics_me_merge_candidates.py",
    "code/src/analyze_me_merge_candidates.py",
    "code/src/run_official_kme_module_diagnosis.py",
    "code/src/run_pairwise_hd_module_diagnosis.py",
    "code/src/run_cross_omics_kme_weighted_jaccard.py",
    "code/src/build_final_jaccard_protein_scope.py",
    "code/src/build_final_jaccard_figure.py",
    "code/src/build_final_jaccard_bubble_figure.py",
    "code/src/final_jaccard_dashboard.html",
    "code/src/run_transcriptomics_outlier_removed_diagnostic.py",
    "code/src/backfill_outlier_after_preprocess_qc.py",
    "code/src/run_transcriptomics_outlier_removed_find_modules.py",
    "code/src/run_transcriptomics_outlier_removed_module_diagnosis.py",
    "code/src/finalize_outlier_removal_no_cut.py",
    "code/src/validate_transcriptomics_outlier_removed_diagnostic.py",
    "code/src/validate_transcriptomics_outlier_removed_downstream.py",
    "code/src/validate_cross_omics_kme_weighted_jaccard.py",
    "code/src/pipeline_provenance.py",
    "code/tests/test_pairwise_hd_module_diagnosis.py",
    "code/tests/test_cross_omics_kme_weighted_jaccard.py",
    "code/tests/test_final_jaccard_bubble_figure.py",
    "code/tests/test_final_jaccard_navigation.py",
    "code/tests/test_final_jaccard_protein_scope.py",
    "code/tests/test_after_preprocess_outlier_navigation.py",
    "code/tests/test_proteomics_qc_comparison.py",
    "PIPELINE_SPEC.md",
    "CHANGELOG.md",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_verified(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    if source.stat().st_size != target.stat().st_size or sha256(source) != sha256(target):
        raise IOError(f"Copy verification failed: {source} -> {target}")


def read_dashboard_bundle(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8").strip()
    if not text.startswith(DATA_PREFIX) or not text.endswith(";"):
        raise ValueError(f"Unexpected dashboard data wrapper: {path}")
    return json.loads(text[len(DATA_PREFIX) : -1])


def safe_part(value: str, fallback: str = "item") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return (cleaned or fallback)[:80]


def _convert_lossless_webp(task: tuple[str, str]) -> dict[str, Any]:
    from PIL import Image

    source = Path(task[0])
    target = Path(task[1])
    Image.MAX_IMAGE_PIXELS = None
    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as image:
        mode = "RGBA" if "A" in image.getbands() else "RGB"
        pixels = image.convert(mode)
        source_pixel_hash = hashlib.sha256(pixels.tobytes()).hexdigest()
        size = pixels.size
        pixels.save(target, "WEBP", lossless=True, method=2, exact=True)
    with Image.open(target) as converted:
        decoded = converted.convert(mode)
        target_pixel_hash = hashlib.sha256(decoded.tobytes()).hexdigest()
        if decoded.size != size or target_pixel_hash != source_pixel_hash:
            raise ValueError(f"Lossless pixel validation failed: {source} -> {target}")
    return {
        "source": source.as_posix(),
        "target": target.as_posix(),
        "source_bytes": source.stat().st_size,
        "target_bytes": target.stat().st_size,
        "width": size[0],
        "height": size[1],
        "pixel_sha256": source_pixel_hash,
    }


def rewrite_plotly_html(source: Path, target: Path) -> int:
    text = source.read_text(encoding="utf-8")
    rewritten, replacements = PLOTLY_CDN_RE.subn(
        '<script charset="utf-8" src="../vendor/plotly.min.js"></script>', text
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(rewritten, encoding="utf-8")
    parser = _RemoteResourceParser()
    parser.feed(rewritten)
    if parser.remote_resources:
        raise ValueError(
            f"Remote browser resource remains in {source}: {parser.remote_resources[:3]}"
        )
    return replacements


class _RemoteResourceParser(HTMLParser):
    """Find resources loaded by actual HTML tags, ignoring URLs inside JS text."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=False)
        self.remote_resources: list[str] = []

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        values = dict(attrs)
        candidates: list[str | None] = []
        if tag in {"script", "img", "iframe", "source", "video", "audio"}:
            candidates.append(values.get("src"))
        if tag == "link":
            candidates.append(values.get("href"))
        for value in candidates:
            if value and value.lower().startswith(("http://", "https://")):
                self.remote_resources.append(value)


def write_readme(path: Path, *, records: int, artifacts: int, files: int) -> None:
    path.write_text(
        f"""# Gene-analysis complete offline inspection dashboard

이 저장소는 `results_20261119`에서 로컬로 보던 experiment inspection dashboard를
다른 컴퓨터에서도 같은 구조로 검토하기 위한 **private offline bundle**입니다.

## 다른 컴퓨터에서 보기

```powershell
git clone https://github.com/sumincho0830/gene-analysis-pairwise-hd-offline.git
cd gene-analysis-pairwise-hd-offline
start index.html
```

별도 Python 환경이나 web server가 필요하지 않습니다. `index.html`을 직접 열면
기본 **Cell type comparison** 화면에서 동일 cell type의 모든 조건을 산출물 탭별로
비교할 수 있습니다. **Individual run inspection**과 review CSV import/export도
로컬 dashboard와 동일하게 유지됩니다.

## 포함된 범위

- experiment records: {records}
- visual/interactive artifacts: {artifacts}
- downloadable CSV/JSON/provenance files: {files}
- feature selection, before/after heatmap·PCA, interactive PCA 3D
- sample-cut diagnostic, soft threshold, gene/ME dendrogram, module sizes
- `Ng_MOG1`·`Ng_HD3` 제거 전/후 결과 비교; 제거 후 branch에는 finite sample cut 없음
- module–diagnosis 및 RNA–protein module similarity
- dashboard와 결과를 재생성하는 핵심 Python 코드와 pipeline 기록

PNG figure는 decoded pixel이 원본과 완전히 같은 lossless WebP로만 재포장했습니다.
Plotly HTML은 저장소에 포함된 `assets/vendor/plotly.min.js`를 사용하므로 인터넷 없이
회전·확대할 수 있습니다.

## 데이터 취급

이 배포본은 현재 로컬 dashboard와 같은 inspection 경험을 보존하므로 interactive
PCA hover와 일부 QC 파일에 source sample ID가 포함될 수 있습니다. 저장소를 private로
유지하고 승인된 사용자에게만 접근 권한을 부여해야 합니다. raw/formatted expression
matrix와 PyWGCNA pickle object는 포함하지 않습니다.

## 무결성 확인

```powershell
python verify_release.py
```

모든 파일이 `SHA256SUMS.txt`와 일치하면 `OK`를 출력합니다.
""",
        encoding="utf-8",
    )


def write_verifier(path: Path) -> None:
    path.write_text(
        """from __future__ import annotations

import hashlib
from pathlib import Path
import sys


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


root = Path(__file__).resolve().parent
failures = []
for line in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    expected, relative = line.split("  ", 1)
    target = root / relative
    if not target.is_file():
        failures.append(f"MISSING {relative}")
    elif sha256(target) != expected:
        failures.append(f"HASH {relative}")

if failures:
    print("FAILED")
    print("\\n".join(failures))
    sys.exit(1)
print("OK")
""",
        encoding="utf-8",
    )


def build(
    run_root: Path,
    output_root: Path,
    *,
    plotly_js: Path,
    workers: int,
    resume: bool,
) -> dict[str, Any]:
    dashboard_root = run_root / "dashboard"
    bundle = read_dashboard_bundle(dashboard_root / "dashboard-data.js")
    if output_root.exists() and any(output_root.iterdir()) and not resume:
        raise FileExistsError(f"Output directory must be empty: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    # This branch was intentionally retired: after removing Ng_MOG1 and Ng_HD3,
    # preprocessing stays at cut=inf.  A resumed release may still contain HTML
    # assets created by the earlier selected-cut prototype, so remove that exact
    # generated directory before computing the release manifest/checksums.
    obsolete_selected_cut = (
        output_root / "assets" / "after_outlier_selected_cut"
    ).resolve()
    if output_root.resolve() not in obsolete_selected_cut.parents:
        raise ValueError(f"Obsolete-asset path escaped release root: {obsolete_selected_cut}")
    if obsolete_selected_cut.is_dir():
        shutil.rmtree(obsolete_selected_cut)

    copy_verified(dashboard_root / "index.html", output_root / "index.html")
    copy_verified(plotly_js, output_root / "assets" / "vendor" / "plotly.min.js")
    publication_files: list[Path] = []
    final_jaccard_page = dashboard_root / "final-jaccard.html"
    if final_jaccard_page.is_file():
        target = output_root / "final-jaccard.html"
        copy_verified(final_jaccard_page, target)
        publication_files.append(target)
    for relative_root in (
        Path("figures/final_jaccard"),
        Path("figures/final_jaccard_bubble"),
    ):
        source_root = dashboard_root / relative_root
        if not source_root.is_dir():
            continue
        for source in source_root.rglob("*"):
            if source.is_file():
                target = output_root / source.relative_to(dashboard_root)
                copy_verified(source, target)
                publication_files.append(target)
    module_member_root = dashboard_root / "module-members"
    module_member_files: list[Path] = []
    if module_member_root.is_dir():
        for source in module_member_root.rglob("*"):
            if source.is_file():
                target = output_root / source.relative_to(dashboard_root)
                copy_verified(source, target)
                module_member_files.append(target)
    cross_overlap_root = dashboard_root / "cross-omics-overlap"
    cross_overlap_files: list[Path] = []
    if cross_overlap_root.is_dir():
        for source in cross_overlap_root.rglob("*"):
            if source.is_file():
                target = output_root / source.relative_to(dashboard_root)
                copy_verified(source, target)
                cross_overlap_files.append(target)

    artifact_sources: dict[str, tuple[Path, Path, str]] = {}
    file_sources: dict[str, tuple[Path, Path]] = {}
    release_records: list[dict[str, Any]] = []

    for record in bundle["records"]:
        release_record = dict(record)
        release_artifacts: list[dict[str, Any]] = []
        for artifact in record.get("artifacts", []):
            source = (dashboard_root / artifact["path"]).resolve()
            if not source.is_file():
                raise FileNotFoundError(source)
            source_key = source.as_posix().lower()
            if source_key not in artifact_sources:
                digest = sha256(source)[:16]
                artifact_id = safe_part(str(artifact.get("id", "artifact")))
                if artifact.get("kind") == "image":
                    target = output_root / "assets" / artifact_id / f"{digest}.webp"
                    operation = "webp"
                elif artifact.get("kind") == "html":
                    target = output_root / "assets" / artifact_id / f"{digest}.html"
                    operation = "html"
                else:
                    target = output_root / "assets" / artifact_id / f"{digest}{source.suffix}"
                    operation = "copy"
                artifact_sources[source_key] = (source, target, operation)
            target = artifact_sources[source_key][1]
            release_artifacts.append(
                {**artifact, "path": target.relative_to(output_root).as_posix()}
            )
        release_record["artifacts"] = release_artifacts

        release_files: list[dict[str, Any]] = []
        for item in record.get("files", []):
            source = (dashboard_root / item["path"]).resolve()
            if not source.is_file():
                raise FileNotFoundError(source)
            source_key = source.as_posix().lower()
            if source_key not in file_sources:
                digest = sha256(source)[:16]
                category = safe_part(str(item.get("category", "data")))
                target = (
                    output_root
                    / "data"
                    / category
                    / f"{digest}_{safe_part(source.name, 'file')}"
                )
                file_sources[source_key] = (source, target)
            target = file_sources[source_key][1]
            release_files.append({**item, "path": target.relative_to(output_root).as_posix()})
        release_record["files"] = release_files
        release_records.append(release_record)

    image_pairs = [
        (source.as_posix(), target.as_posix())
        for source, target, operation in artifact_sources.values()
        if operation == "webp"
    ]
    image_tasks = [task for task in image_pairs if not Path(task[1]).is_file()]
    image_results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = [executor.submit(_convert_lossless_webp, task) for task in image_tasks]
        for index, future in enumerate(as_completed(futures), start=1):
            image_results.append(future.result())
            if index % 100 == 0 or index == len(futures):
                print(f"Converted {index}/{len(futures)} lossless figures", flush=True)

    plotly_replacements = 0
    for source, target, operation in artifact_sources.values():
        if operation == "html":
            plotly_replacements += rewrite_plotly_html(source, target)
        elif operation == "copy":
            copy_verified(source, target)
    for source, target in file_sources.values():
        copy_verified(source, target)

    release_bundle = {
        **bundle,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": "offline_full_dashboard_release",
        "sources": [
            {**source, "summary_path": "", "record_count": sum(
                record["source_id"] == source["id"] for record in release_records
            )}
            for source in bundle["sources"]
        ],
        "records": release_records,
    }
    (output_root / "dashboard-data.js").write_text(
        DATA_PREFIX
        + json.dumps(release_bundle, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )

    for relative in CODE_FILES:
        source = run_root / relative
        if source.is_file():
            copy_verified(source, output_root / "reproducibility" / relative)
    write_readme(
        output_root / "README.md",
        records=len(release_records),
        artifacts=sum(len(record.get("artifacts", [])) for record in release_records),
        files=sum(len(record.get("files", [])) for record in release_records),
    )
    write_verifier(output_root / "verify_release.py")
    (output_root / ".gitattributes").write_text("* -text\n", encoding="utf-8")
    (output_root / ".gitignore").write_text(
        "__pycache__/\n*.pyc\n", encoding="utf-8"
    )

    source_image_bytes = sum(Path(source).stat().st_size for source, _ in image_pairs)
    unique_image_targets = {Path(target) for _, target in image_pairs}
    release_image_bytes = sum(target.stat().st_size for target in unique_image_targets)
    manifest = {
        "schema_version": 2,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": "results_20261119",
        "records": len(release_records),
        "artifact_references": sum(
            len(record.get("artifacts", [])) for record in release_records
        ),
        "unique_artifacts": len(artifact_sources),
        "file_references": sum(len(record.get("files", [])) for record in release_records),
        "unique_files": len(file_sources),
        "module_member_files": len(module_member_files),
        "cross_omics_overlap_files": len(cross_overlap_files),
        "publication_files": len(publication_files),
        "lossless_webp_sources": len(image_pairs),
        "lossless_webp_files": len(unique_image_targets),
        "lossless_webp_converted_this_run": len(image_results),
        "source_image_bytes": source_image_bytes,
        "release_image_bytes": release_image_bytes,
        "image_size_reduction_fraction": (
            1 - release_image_bytes / source_image_bytes if source_image_bytes else 0
        ),
        "plotly_cdn_replacements": plotly_replacements,
        "remote_browser_resources": 0,
        "privacy_scope": {
            "raw_expression_matrices": "excluded",
            "pywgcna_objects": "excluded",
            "dashboard_linked_outputs": "included",
            "source_sample_ids_in_interactive_qc": "included_private_repo_only",
        },
    }
    (output_root / "release_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    hash_paths = sorted(
        path
        for path in output_root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    (output_root / "SHA256SUMS.txt").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(output_root).as_posix()}"
            for path in hash_paths
        )
        + "\n",
        encoding="utf-8",
    )
    manifest["release_files"] = len(hash_paths) + 1
    (output_root / "release_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    # Manifest changed after the first checksum pass; rebuild the checksums once.
    hash_paths = sorted(
        path
        for path in output_root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    (output_root / "SHA256SUMS.txt").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(output_root).as_posix()}"
            for path in hash_paths
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    default_run_root = Path(__file__).resolve().parents[2]
    parser.add_argument("--run-root", type=Path, default=default_run_root)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=default_run_root.parent / "results_20261119_full_dashboard_staging",
    )
    parser.add_argument(
        "--plotly-js",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / ".venv"
        / "Lib"
        / "site-packages"
        / "plotly"
        / "package_data"
        / "plotly.min.js",
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = build(
        args.run_root.resolve(),
        args.output_root.resolve(),
        plotly_js=args.plotly_js.resolve(),
        workers=args.workers,
        resume=args.resume,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
