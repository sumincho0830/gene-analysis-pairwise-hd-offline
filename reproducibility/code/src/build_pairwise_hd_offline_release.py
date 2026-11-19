"""Build a privacy-minimized offline release of pairwise-HD results.

The release contains the static dashboard, pairwise disease-vs-HD heatmaps,
aggregate module-level tables, manifests, and the relevant source code. Raw
expression matrices, individual sample identifiers, and PyWGCNA pickle objects
are intentionally excluded.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any


DATA_PREFIX = "window.EXPERIMENT_DASHBOARD_DATA = "
PAIRWISE_FILE_IDS = {
    "diagnosis_pairwise_statistics",
    "diagnosis_pairwise_correlations",
    "diagnosis_pairwise_pvalues",
    "diagnosis_pairwise_qvalues",
    "diagnosis_pairwise_counts",
    "diagnosis_pairwise_summary",
    "diagnosis_pairwise_manifest",
}
CODE_FILES = (
    "code/src/run_pairwise_hd_module_diagnosis.py",
    "code/src/run_official_kme_module_diagnosis.py",
    "code/src/pipeline_provenance.py",
    "code/src/build_experiment_inspection_dashboard.py",
    "code/src/build_pairwise_hd_offline_release.py",
    "code/tests/test_pairwise_hd_module_diagnosis.py",
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


def safe_record_key(record_id: str) -> str:
    prefix = hashlib.sha256(record_id.encode("utf-8")).hexdigest()[:12]
    readable = "_".join(
        part.replace(" ", "_")[:24] for part in record_id.split("::")[-3:]
    )
    return f"{prefix}_{readable}"[:90]


def sanitize_manifest(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if "run_root" in payload:
        payload["run_root"] = "source run: results_20261119 (not included in release)"
    return payload


def write_readme(path: Path, record_count: int) -> None:
    path.write_text(
        f"""# Pairwise disease-vs-HD module diagnosis offline release

이 저장소는 `results_20261119`의 **집계된 module–diagnosis 결과만** 오프라인에서 검토하기 위한 배포본입니다.

## 다른 컴퓨터에서 보기

```powershell
git clone <PRIVATE_REPOSITORY_URL>
cd <CLONED_DIRECTORY>
start index.html
```

별도 Python 환경이나 web server가 필요하지 않습니다. 인터넷 연결을 끊은 상태에서도 `index.html`을 직접 열 수 있습니다.

## 포함된 결과

- 실험 record: {record_count}개
- `MS vs HD`, `MOGAD vs HD`, `NMOSD vs HD` heatmap
- signed point-biserial `r`, two-sided `p`, within-network BH `q`
- 질환별 집계 sample 수, module feature 수, 분석 summary와 provenance
- 결과를 생성한 Python 코드와 단위 테스트

## 의도적으로 제외한 데이터

- raw/formatted expression matrix
- 개인 sample ID가 포함된 sample-level table
- PyWGCNA `.p`, `.p.gz`, `.p.zst` 객체
- 전체 feature-level kME table

따라서 이 저장소는 **결과 열람과 방법 감사용**입니다. 분석을 처음부터 재실행하려면 별도로 보관된 원본 데이터와 전체 `results_20261119` 실행 디렉토리가 필요합니다.

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


root = Path(__file__).resolve().parent
failures = []
for line in (root / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
    expected, relative = line.split("  ", 1)
    target = root / relative
    if not target.is_file():
        failures.append(f"MISSING {relative}")
        continue
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    if digest != expected:
        failures.append(f"HASH {relative}")

if failures:
    print("FAILED")
    print("\\n".join(failures))
    sys.exit(1)
print("OK")
""",
        encoding="utf-8",
    )


def build(run_root: Path, output_root: Path) -> dict[str, Any]:
    dashboard_root = run_root / "dashboard"
    bundle = read_dashboard_bundle(dashboard_root / "dashboard-data.js")
    output_root.mkdir(parents=True, exist_ok=True)

    copy_verified(dashboard_root / "index.html", output_root / "index.html")
    release_records: list[dict[str, Any]] = []
    copied_outputs = 0

    for record in bundle["records"]:
        release_record = dict(record)
        pairwise_artifacts = [
            artifact
            for artifact in record.get("artifacts", [])
            if artifact.get("id") == "diagnosis_pairwise_hd"
        ]
        if not pairwise_artifacts:
            # Three dashboard records did not complete findModules and therefore
            # have no module eigengenes or pairwise diagnosis result to release.
            continue
        if len(pairwise_artifacts) != 1:
            raise ValueError(f"Expected one pairwise artifact: {record['id']}")
        source_artifact = (dashboard_root / pairwise_artifacts[0]["path"]).resolve()
        artifact_target = output_root / "assets" / "pairwise_hd" / source_artifact.name
        copy_verified(source_artifact, artifact_target)
        release_record["artifacts"] = [
            {
                **pairwise_artifacts[0],
                "path": artifact_target.relative_to(output_root).as_posix(),
            }
        ]

        record_key = safe_record_key(str(record["id"]))
        release_files = []
        for item in record.get("files", []):
            if item.get("id") not in PAIRWISE_FILE_IDS:
                continue
            source = (dashboard_root / item["path"]).resolve()
            target = output_root / "data" / record_key / source.name
            if item.get("id") == "diagnosis_pairwise_manifest":
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(
                    json.dumps(sanitize_manifest(source), indent=2, ensure_ascii=False),
                    encoding="utf-8",
                )
            else:
                copy_verified(source, target)
            release_files.append({**item, "path": target.relative_to(output_root).as_posix()})
            copied_outputs += 1
        if len(release_files) != len(PAIRWISE_FILE_IDS):
            raise ValueError(
                f"Pairwise file set is incomplete for {record['id']}: {len(release_files)}"
            )
        release_record["files"] = release_files
        release_records.append(release_record)

    release_bundle = {
        **bundle,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": "offline_pairwise_hd_release",
        "sources": [
            {**source, "summary_path": "", "record_count": sum(
                record["source_id"] == source["id"] for record in release_records
            )}
            for source in bundle["sources"]
        ],
        "records": release_records,
    }
    (output_root / "dashboard-data.js").write_text(
        DATA_PREFIX + json.dumps(release_bundle, ensure_ascii=False, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )

    for relative in CODE_FILES:
        copy_verified(run_root / relative, output_root / "reproducibility" / relative)
    write_readme(output_root / "README.md", len(release_records))
    write_verifier(output_root / "verify_release.py")

    manifest = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_run": "results_20261119",
        "records": len(release_records),
        "pairwise_heatmaps": len(release_records),
        "pairwise_output_files": copied_outputs,
        "privacy_scope": {
            "sample_level_identifiers": "excluded",
            "expression_matrices": "excluded",
            "pywgcna_objects": "excluded",
            "feature_level_kme": "excluded",
            "module_level_aggregate_statistics": "included",
        },
    }
    (output_root / "release_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    hash_paths = sorted(
        path
        for path in output_root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    (output_root / "SHA256SUMS.txt").write_text(
        "\n".join(
            f"{sha256(path)}  {path.relative_to(output_root).as_posix()}" for path in hash_paths
        )
        + "\n",
        encoding="utf-8",
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(__file__).resolve().parents[3] / "results_20261119_pairwise_hd_offline",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = build(args.run_root.resolve(), args.output_root.resolve())
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
