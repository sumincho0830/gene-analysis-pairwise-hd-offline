"""Build an offline, extensible experiment-inspection dashboard.

The dashboard UI is stable. Experiment-specific adapters normalize saved
results into a shared record/artifact/provenance schema and emit a local
``dashboard-data.js`` bundle, so direct ``file://`` use does not require fetch
or a web server.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any


SCHEMA_VERSION = "1.0"
MODULE_CONDITION = "01_default_ds2_mms50_med0p20"


# Visual layer derived from results_20261015/design-system.md.  The dashboard
# keeps its existing data schema and interactions; this stylesheet only
# changes hierarchy, spacing, typography, and responsive behavior.
DESIGN_SYSTEM_CSS = r"""
/* Scientific report design system --------------------------------------- */
:root {
  --font-ui: Arial, Helvetica, "Malgun Gothic", sans-serif;
  --surface: #ffffff;
  --surface-subtle: #fafafa;
  --text: #262626;
  --text-secondary: #555555;
  --text-muted: #6a6a6a;
  --link: #353535;
  --line: #dddddd;
  --line-strong: #b8b8b8;
  --line-subtle: #eeeeee;
  --hover: #f5f5f5;
  --selected: #f1f1f1;
  --ok: #2f6b4f;
  --warn: #955b2a;
  --page-max: 1560px;
  --list-width: 360px;
}

* { box-sizing: border-box; }
html { color-scheme: light; background: var(--surface); }
body {
  margin: 0;
  background: var(--surface);
  color: var(--text);
  font-family: var(--font-ui);
  font-size: 16px;
  font-weight: 400;
  line-height: 1.55;
  -webkit-font-smoothing: antialiased;
  text-rendering: optimizeLegibility;
}

button, input, select, textarea { font: inherit; color: inherit; }
button, select, input, textarea { border-radius: 0; }
button:focus-visible, a:focus-visible, select:focus-visible,
input:focus-visible, textarea:focus-visible {
  outline: 2px solid #555;
  outline-offset: 2px;
}

.top {
  max-width: var(--page-max);
  margin: 0 auto;
  padding: 38px 32px 24px;
  border-bottom: 1px solid var(--line-subtle);
  text-align: center;
  background: var(--surface);
}
.top h1 {
  margin: 0 0 7px;
  font-size: 1.55rem;
  line-height: 1.4;
  font-weight: 400;
  letter-spacing: 0;
}
.top .sub {
  margin: 0;
  color: var(--text-muted);
  font-size: .82rem;
  line-height: 1.6;
}

.kpis {
  max-width: var(--page-max);
  margin: 0 auto;
  display: grid;
  grid-template-columns: repeat(6, minmax(100px, 1fr));
  gap: 0;
  background: var(--surface);
  border-bottom: 1px solid var(--line);
}
.kpi {
  min-width: 0;
  padding: 17px 20px 15px;
  background: var(--surface);
  border-right: 1px solid var(--line-subtle);
  text-align: center;
}
.kpi:last-child { border-right: 0; }
.kpi b {
  display: block;
  margin-bottom: 1px;
  color: var(--text);
  font-size: 1.25rem;
  line-height: 1.25;
  font-weight: 400;
  font-variant-numeric: tabular-nums;
}
.kpi span {
  color: var(--text-muted);
  font-size: .78rem;
  line-height: 1.4;
}

.toolbar {
  position: sticky;
  top: 0;
  z-index: 20;
  max-width: var(--page-max);
  margin: 0 auto;
  padding: 11px 16px;
  display: flex;
  align-items: center;
  gap: 7px;
  flex-wrap: wrap;
  background: rgba(255,255,255,.97);
  border-bottom: 1px solid var(--line);
  backdrop-filter: blur(8px);
}
.toolbar select, .toolbar input {
  min-width: 125px;
  height: 34px;
  padding: 6px 30px 6px 9px;
  border: 1px solid var(--line);
  background-color: var(--surface);
  color: var(--text-secondary);
  font-size: .87rem;
  line-height: 1.2;
}
.toolbar select {
  appearance: none;
  background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='8' viewBox='0 0 12 8'%3E%3Cpath d='M1 1.5 6 6.5l5-5' fill='none' stroke='%23666' stroke-width='1.25'/%3E%3C/svg%3E");
  background-repeat: no-repeat;
  background-position: right 9px center;
}
.toolbar #search {
  min-width: 215px;
  padding-right: 9px;
}
.toolbar input::placeholder { color: #8a8a8a; }
.toolbar button, .detail-head button, .compare button {
  min-height: 34px;
  padding: 6px 10px;
  border: 1px solid var(--line-strong);
  background: var(--surface);
  color: var(--link);
  font-size: .84rem;
  line-height: 1.2;
  cursor: pointer;
}
.toolbar button:hover, .detail-head button:hover, .compare button:hover {
  background: var(--hover);
}
.spacer { flex: 1 1 12px; }

.layout {
  max-width: var(--page-max);
  margin: 0 auto;
  display: grid;
  grid-template-columns: var(--list-width) minmax(0, 1fr);
  min-height: calc(100vh - 235px);
  border-left: 1px solid var(--line-subtle);
  border-right: 1px solid var(--line-subtle);
}
.list {
  height: calc(100vh - 235px);
  overflow: auto;
  background: var(--surface);
  border-right: 1px solid var(--line);
  scrollbar-width: thin;
  scrollbar-color: #c7c7c7 transparent;
}
.list-head {
  position: sticky;
  top: 0;
  z-index: 2;
  padding: 11px 16px 10px;
  background: rgba(255,255,255,.97);
  border-bottom: 1px solid var(--line-strong);
  color: var(--text-muted);
  font-size: .78rem;
  font-variant-numeric: tabular-nums;
}
.record {
  position: relative;
  width: 100%;
  padding: 13px 16px 12px 18px;
  text-align: left;
  background: var(--surface);
  border: 0;
  border-bottom: 1px solid var(--line-subtle);
  cursor: pointer;
}
.record::before {
  content: "";
  position: absolute;
  inset: 0 auto 0 0;
  width: 2px;
  background: transparent;
}
.record:hover { background: var(--hover); }
.record.active { background: var(--selected); }
.record.active::before { background: var(--text); }
.record-title {
  color: var(--text);
  font-size: .87rem;
  line-height: 1.4;
  font-weight: 500;
}
.record-meta, .record-metrics {
  display: flex;
  gap: 5px 9px;
  flex-wrap: wrap;
  margin-top: 4px;
  color: var(--text-muted);
  font-size: .77rem;
  line-height: 1.45;
  font-variant-numeric: tabular-nums;
}
.pill {
  display: inline;
  padding: 0;
  border: 0;
  border-radius: 0;
  color: var(--text-muted);
}
.pill.warn { color: var(--warn); }
.pill.ok { color: var(--ok); }

.detail {
  min-width: 0;
  padding: 30px 36px 56px;
  background: var(--surface);
}
.detail-head {
  display: flex;
  align-items: flex-start;
  gap: 20px;
  padding-bottom: 15px;
  border-bottom: 1px solid var(--line-strong);
}
.detail-head h2 {
  margin: 0 0 4px;
  font-size: 1.32rem;
  line-height: 1.4;
  font-weight: 400;
}
.detail-head .sub {
  color: var(--text-muted);
  font-size: .81rem;
  line-height: 1.5;
}
.detail-head .actions { margin-left: auto; }

.flags {
  display: flex;
  gap: 7px 18px;
  flex-wrap: wrap;
  margin: 11px 0 15px;
}
.flag {
  padding: 0;
  border: 0;
  background: transparent;
  color: var(--warn);
  font-size: .81rem;
  line-height: 1.5;
}
.flag::before { content: "Status · "; color: var(--text-muted); }

.metric-grid {
  display: grid;
  grid-template-columns: repeat(6, minmax(110px, 1fr));
  gap: 0;
  margin: 0 0 28px;
  background: var(--surface);
  border-top: 1px solid var(--line-subtle);
  border-bottom: 1px solid var(--line-subtle);
}
.metric {
  min-width: 0;
  padding: 11px 12px 10px 0;
  background: var(--surface);
  border-right: 1px solid var(--line-subtle);
}
.metric:nth-child(6n) { border-right: 0; }
.metric .label {
  margin-bottom: 2px;
  color: var(--text-muted);
  font-size: .78rem;
  line-height: 1.35;
}
.metric .value {
  color: var(--text);
  font-size: .96rem;
  line-height: 1.35;
  font-weight: 400;
  font-variant-numeric: tabular-nums;
  overflow-wrap: anywhere;
}

.review {
  display: grid;
  grid-template-columns: minmax(290px, .85fr) minmax(340px, 1.15fr);
  gap: 36px;
  margin: 0 0 34px;
  padding: 18px 0 20px;
  background: var(--surface);
  border: 0;
  border-top: 1px solid var(--line-strong);
  border-bottom: 1px solid var(--line-subtle);
  font-size: .88rem;
}
.review b {
  font-size: .88rem;
  font-weight: 500;
}
.decision-row {
  display: flex;
  gap: 20px;
  margin: 9px 0 12px;
  border-bottom: 1px solid var(--line-subtle);
}
.decision-row button {
  margin: 0 0 -1px;
  padding: 5px 0 7px;
  border: 0;
  border-bottom: 2px solid transparent;
  background: transparent;
  color: var(--text-muted);
  cursor: pointer;
}
.decision-row button:hover { color: var(--text); }
.decision-row button.selected {
  border-color: var(--text);
  background: transparent;
  color: var(--text);
}
.checklist label {
  display: block;
  margin: 6px 0;
  color: var(--text-secondary);
  font-size: .81rem;
}
.checklist input { accent-color: #333; }
.review textarea {
  width: 100%;
  min-height: 96px;
  margin-top: 9px;
  padding: 9px 10px;
  resize: vertical;
  border: 1px solid var(--line);
  background: var(--surface);
  color: var(--text);
  font-size: .84rem;
  line-height: 1.55;
}
.review .sub { margin-top: 5px; font-size: .76rem; color: var(--text-muted); }

.section-title {
  margin: 28px 0 11px;
  color: var(--text);
  font-size: .92rem;
  line-height: 1.45;
  font-weight: 500;
}
.artifact-tabs {
  display: flex;
  align-items: flex-end;
  gap: 22px;
  overflow-x: auto;
  padding: 0;
  border-bottom: 1px solid var(--line-subtle);
  scrollbar-width: thin;
}
.artifact-tabs button {
  flex: 0 0 auto;
  margin: 0 0 -1px;
  padding: 7px 0 9px;
  border: 0;
  border-bottom: 2px solid transparent;
  background: transparent;
  color: var(--text-muted);
  font-size: .82rem;
  line-height: 1.3;
  cursor: pointer;
  white-space: nowrap;
}
.artifact-tabs button:hover { color: var(--text); }
.artifact-tabs button.active {
  border-color: var(--text);
  background: transparent;
  color: var(--text);
}
.viewer {
  min-height: 460px;
  padding: 22px 0 10px;
  display: flex;
  align-items: flex-start;
  justify-content: center;
  overflow: auto;
  background: var(--surface);
  border: 0;
  border-bottom: 1px solid var(--line-subtle);
}
.viewer > a { display: block; width: 100%; text-align: center; }
.viewer img {
  display: block;
  max-width: 100%;
  height: auto;
  margin: 0 auto;
}
.viewer iframe { width: 100%; height: 720px; border: 1px solid var(--line-subtle); }
.viewer .open { padding: 0 0 12px; font-size: .82rem; }
.viewer a, .files a {
  color: var(--link);
  text-decoration: underline;
  text-decoration-color: #aaa;
  text-underline-offset: 3px;
}

.lower {
  display: grid;
  grid-template-columns: minmax(0, 1.1fr) minmax(320px, .9fr);
  gap: 44px;
  align-items: start;
}
.table-wrap {
  max-height: 420px;
  overflow: auto;
  border: 0;
  border-top: 1px solid var(--line-strong);
}
table { width: 100%; border-collapse: collapse; font-size: .82rem; line-height: 1.55; }
th, td {
  padding: 8px 10px 8px 0;
  border-bottom: 1px solid var(--line-subtle);
  text-align: left;
  white-space: nowrap;
  font-weight: 400;
  font-variant-numeric: tabular-nums;
}
th {
  position: sticky;
  top: 0;
  z-index: 1;
  background: var(--surface);
  color: var(--text-muted);
  font-size: .77rem;
  border-bottom-color: var(--line-strong);
}
.swatch {
  display: inline-block;
  width: 10px;
  height: 10px;
  margin-right: 7px;
  vertical-align: -1px;
  border: 1px solid #777;
}
.files {
  padding: 0;
  border: 0;
  border-top: 1px solid var(--line-strong);
  font-size: .81rem;
}
.files details { padding: 8px 0; border-bottom: 1px solid var(--line-subtle); }
.files summary { color: var(--text-secondary); cursor: pointer; }
.files a { display: block; padding: 4px 0 1px; word-break: break-all; }

.compare {
  max-width: var(--page-max);
  margin: 0 auto;
  padding: 24px 32px 30px;
  border: 0;
  border-bottom: 1px solid var(--line);
  background: var(--surface);
}
.compare-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  gap: 28px;
}
.compare-card { padding: 0; border: 0; border-top: 1px solid var(--line-strong); }
.compare-card h3 { margin: 10px 0 4px; font-size: .84rem; font-weight: 500; }
.compare-card img { width: 100%; height: auto; margin-top: 10px; }
.empty { padding: 70px 20px; color: var(--text-muted); text-align: center; }

/* Cell-centric comparison ------------------------------------------------ */
body[data-view="cell"] .run-only { display: none !important; }
body[data-view="run"] .cell-only { display: none !important; }
.cell-comparison {
  max-width: var(--page-max);
  margin: 0 auto;
  padding: 31px 32px 64px;
  border-left: 1px solid var(--line-subtle);
  border-right: 1px solid var(--line-subtle);
}
.comparison-head {
  display: flex;
  align-items: flex-end;
  justify-content: space-between;
  gap: 24px;
  margin-bottom: 18px;
}
.comparison-head h2 {
  margin: 0 0 4px;
  font-size: 1.35rem;
  line-height: 1.35;
  font-weight: 400;
}
.comparison-counts {
  color: var(--text-muted);
  font-size: .8rem;
  font-variant-numeric: tabular-nums;
  text-align: right;
}
.comparison-note {
  margin: 11px 0 24px;
  color: var(--text-muted);
  font-size: .79rem;
}
.cell-comparison .artifact-tabs {
  flex-wrap: wrap;
  gap: 3px 20px;
  overflow: visible;
}
.condition-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 34px 28px;
  align-items: start;
}
.condition-panel {
  min-width: 0;
  padding-top: 10px;
  border-top: 1px solid var(--line-strong);
}
.condition-head {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: 14px;
  min-height: 56px;
  margin-bottom: 9px;
}
.condition-head h3 {
  margin: 0 0 3px;
  font-size: .91rem;
  line-height: 1.4;
  font-weight: 500;
}
.condition-meta {
  color: var(--text-muted);
  font-size: .76rem;
  line-height: 1.5;
  font-variant-numeric: tabular-nums;
}
.condition-review {
  flex: 0 0 auto;
  color: var(--text-muted);
  font-size: .72rem;
  text-transform: uppercase;
}
.condition-review.use { color: var(--ok); }
.condition-review.hold, .condition-review.reject { color: var(--warn); }
.condition-viewer {
  height: 370px;
  display: flex;
  align-items: center;
  justify-content: center;
  overflow: auto;
  background: var(--surface-subtle);
  border: 1px solid var(--line-subtle);
}
.condition-viewer > a { width: 100%; height: 100%; display: flex; align-items: center; justify-content: center; }
.condition-viewer img { display: block; max-width: 100%; max-height: 100%; width: auto; height: auto; }
.condition-viewer iframe { width: 100%; height: 100%; border: 0; background: var(--surface); }
.condition-missing { color: var(--text-muted); font-size: .8rem; text-align: center; padding: 30px; }
.condition-links {
  display: flex;
  gap: 16px;
  align-items: center;
  margin-top: 8px;
  font-size: .77rem;
}
.condition-links a, .condition-links button {
  padding: 0 0 2px;
  border: 0;
  border-bottom: 1px solid #aaa;
  background: transparent;
  color: var(--link);
  cursor: pointer;
  text-decoration: none;
}
.condition-links a:hover, .condition-links button:hover { border-bottom-color: var(--text); }

@media (max-width: 1180px) {
  :root { --list-width: 320px; }
  .metric-grid { grid-template-columns: repeat(4, minmax(105px, 1fr)); }
  .metric:nth-child(6n) { border-right: 1px solid var(--line-subtle); }
  .metric:nth-child(4n) { border-right: 0; }
  .lower { grid-template-columns: 1fr; gap: 24px; }
  .condition-grid { grid-template-columns: 1fr; }
}

@media (max-width: 820px) {
  .top { padding: 26px 16px 19px; }
  .top h1 { font-size: 1.3rem; }
  .kpis { grid-template-columns: repeat(3, 1fr); }
  .kpi:nth-child(3) { border-right: 0; }
  .kpi:nth-child(-n+3) { border-bottom: 1px solid var(--line-subtle); }
  .toolbar { position: static; padding: 10px 12px; }
  .toolbar select, .toolbar input, .toolbar #search { min-width: calc(50% - 4px); flex: 1 1 calc(50% - 4px); }
  .toolbar button { flex: 1 1 auto; }
  .spacer { display: none; }
  .layout { grid-template-columns: 1fr; border-left: 0; border-right: 0; }
  .list { height: 300px; border-right: 0; border-bottom: 1px solid var(--line); }
  .detail { padding: 25px 18px 46px; }
  .review { grid-template-columns: 1fr; gap: 22px; }
  .metric-grid { grid-template-columns: repeat(3, minmax(90px, 1fr)); }
  .metric:nth-child(4n) { border-right: 1px solid var(--line-subtle); }
  .metric:nth-child(3n) { border-right: 0; }
  .viewer { min-height: 300px; }
  .cell-comparison { padding: 25px 18px 48px; border-left: 0; border-right: 0; }
  .comparison-head { display: block; }
  .comparison-counts { margin-top: 8px; text-align: left; }
  .condition-viewer { height: 460px; }
}

@media (max-width: 520px) {
  body { font-size: 15px; }
  .toolbar select, .toolbar input, .toolbar #search { min-width: 100%; flex-basis: 100%; }
  .detail-head { display: block; }
  .detail-head .actions { margin: 12px 0 0; }
  .metric-grid { grid-template-columns: repeat(2, minmax(90px, 1fr)); }
  .metric:nth-child(3n) { border-right: 1px solid var(--line-subtle); }
  .metric:nth-child(2n) { border-right: 0; }
  .artifact-tabs { gap: 17px; }
  .condition-viewer { height: 340px; }
}

@media print {
  .toolbar, .list, .review, .detail-head .actions, .compare { display: none !important; }
  .top { padding: 0 0 20px; }
  .kpis { break-inside: avoid; }
  .layout { display: block; border: 0; }
  .detail { padding: 20px 0 0; }
  .viewer { min-height: 0; }
}
"""


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_number(value: Any) -> int | float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number.is_integer() else number


def relative_path(path: Path, dashboard_root: Path) -> str:
    return Path(path).resolve().relative_to(dashboard_root.parent.resolve()).as_posix()


def dashboard_relative(path: Path, dashboard_root: Path) -> str:
    return "../" + relative_path(path, dashboard_root)


def existing_artifact(
    artifacts: list[dict[str, str]],
    root: Path,
    dashboard_root: Path,
    artifact_id: str,
    label: str,
    relative: str,
    stage: str,
    kind: str = "image",
) -> None:
    path = root / relative
    if path.exists():
        artifacts.append(
            {
                "id": artifact_id,
                "label": label,
                "stage": stage,
                "kind": kind,
                "path": dashboard_relative(path, dashboard_root),
            }
        )


def existing_file(
    files: list[dict[str, str]],
    root: Path,
    dashboard_root: Path,
    file_id: str,
    label: str,
    relative: str,
    category: str,
) -> None:
    path = root / relative
    if path.exists():
        files.append(
            {
                "id": file_id,
                "label": label,
                "category": category,
                "path": dashboard_relative(path, dashboard_root),
            }
        )


def mirrored_image_artifact(
    artifacts: list[dict[str, str]],
    source: Path,
    dashboard_root: Path,
    artifact_id: str,
    label: str,
    stage: str,
) -> None:
    """Mirror long-path images under the dashboard for reliable file:// loading."""

    if not source.exists():
        return
    digest = hashlib.sha256(source.resolve().as_posix().encode("utf-8")).hexdigest()[:20]
    target = dashboard_root / "assets" / "pairwise_hd" / f"{digest}{source.suffix.lower()}"
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or sha256(source) != sha256(target):
        shutil.copy2(source, target)
    artifacts.append(
        {
            "id": artifact_id,
            "label": label,
            "stage": stage,
            "kind": "image",
            "path": dashboard_relative(target, dashboard_root),
        }
    )


def build_variance_filtered_transcriptomics(
    run_root: Path,
    dashboard_root: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    stage_root = run_root / "experiments/04_variance_filtered/transcriptomics"
    summary_path = stage_root / "find_modules_summary.csv"
    rows = read_csv_rows(summary_path)
    records: list[dict[str, Any]] = []

    for row in rows:
        analysis_path = stage_root / row["summary"]
        module_root = analysis_path.parent
        network_root = module_root.parent.parent
        analysis = read_json(analysis_path)
        preprocess_path = network_root / "03_preprocess/02_selected_cut/preprocess_summary.json"
        preprocess = read_json(preprocess_path) if preprocess_path.exists() else {}
        diagnosis_root = network_root / "06_kme_module_diagnosis" / MODULE_CONDITION
        diagnosis_summary_path = diagnosis_root / "analysis_summary.json"
        diagnosis_summary = read_json(diagnosis_summary_path) if diagnosis_summary_path.exists() else {}
        pairwise_summary_path = diagnosis_root / "pairwise_hd_analysis_summary.json"
        pairwise_summary = read_json(pairwise_summary_path) if pairwise_summary_path.exists() else {}
        cross_network_root = (
            run_root
            / "experiments/05_cross_omics_module_similarity/transcriptomics"
            / row["threshold"]
            / row["strategy"]
            / row["mode"]
            / row["specimen"]
        )
        if row["mode"] == "cell_separate":
            cross_network_root = cross_network_root / row["cell_type"]
        cross_summary_path = cross_network_root / "analysis_summary.json"
        cross_summary = read_json(cross_summary_path) if cross_summary_path.exists() else {}
        diagnosis_rows = read_csv_rows(
            diagnosis_root / "module_diagnosis_pairwise_hd_statistics.csv"
        )
        diagnosis_q_values = [
            value
            for value in (
                parse_number(item.get("bh_fdr_q_within_network")) for item in diagnosis_rows
            )
            if value is not None
        ]
        non_estimable_diagnosis = sum(
            parse_number(item.get("p_value_two_sided")) is None for item in diagnosis_rows
        )

        modules = []
        for module in read_csv_rows(module_root / "module_summary.csv"):
            modules.append(
                {
                    "color": module.get("module_color"),
                    "features": parse_number(module.get("feature_count")),
                    "assigned_features": parse_number(module.get("assigned_feature_count")),
                    "unassigned_features": parse_number(module.get("raw_unassigned_feature_count")),
                    "median_abs_kme": parse_number(module.get("assigned_abs_kME_median")),
                    "fully_unassigned": str(module.get("is_fully_raw_unassigned", "")).lower()
                    == "true",
                }
            )

        selected_power = parse_number(analysis.get("selected_power"))
        selected_sft = None
        for metric in read_csv_rows(module_root / "soft_threshold_table.csv"):
            if parse_number(metric.get("Power")) == selected_power:
                selected_sft = {
                    "r2": parse_number(metric.get("SFT.R.sq")),
                    "slope": parse_number(metric.get("slope")),
                    "mean_connectivity": parse_number(metric.get("mean(k)")),
                }
                break

        artifacts: list[dict[str, str]] = []
        mirrored_image_artifact(
            artifacts,
            diagnosis_root / "figures/module_diagnosis_pairwise_hd_heatmap.png",
            dashboard_root,
            "diagnosis_pairwise_hd",
            "Module-diagnosis (vs HD)",
            "module diagnosis",
        )
        for spec in (
            ("me_dendrogram", "ME dendrogram", "05_find_modules/" + MODULE_CONDITION + "/figures/module_eigengene_dendrogram.png", "findModules", "image"),
            ("gene_dendrogram", "Gene dendrogram", "05_find_modules/" + MODULE_CONDITION + "/figures/gene_dendrogram_modules.png", "findModules", "image"),
            ("soft_threshold", "Soft threshold", "05_find_modules/" + MODULE_CONDITION + "/figures/soft_threshold.png", "findModules", "image"),
            ("module_sizes", "Module sizes", "05_find_modules/" + MODULE_CONDITION + "/figures/module_sizes.png", "findModules", "image"),
            ("diagnosis_heatmap", "One-hot diagnosis (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/figures/module_diagnosis_heatmap.png", "module diagnosis", "image"),
            ("feature_selection", "Feature selection", "01_feature_selection/feature_selection.png", "feature selection", "image"),
            ("before_heatmap", "Before heatmap", "02_before_preprocess/02_data_distribution/heatmap.png", "before preprocess", "image"),
            ("before_pca", "Before PCA", "02_before_preprocess/02_data_distribution/pca_metadata.png", "before preprocess", "image"),
            ("before_pca_3d", "Before PCA 3D", "02_before_preprocess/02_data_distribution/pca_3d_interactive.html", "before preprocess", "html"),
            ("cut_diagnostic", "Cut diagnostic", "03_preprocess/01_cut_diagnostic/dendrogram_cut_diagnostic.png", "preprocess", "image"),
            ("selected_cut", "Selected cut", "03_preprocess/02_selected_cut/dendrogram_selected_cut.png", "preprocess", "image"),
            ("after_heatmap", "After heatmap", "04_after_preprocess/02_data_distribution/heatmap.png", "after preprocess", "image"),
            ("after_pca", "After PCA", "04_after_preprocess/02_data_distribution/pca_metadata.png", "after preprocess", "image"),
            ("after_pca_3d", "After PCA 3D", "04_after_preprocess/02_data_distribution/pca_3d_interactive.html", "after preprocess", "html"),
        ):
            existing_artifact(artifacts, network_root, dashboard_root, *spec)
        existing_artifact(
            artifacts,
            cross_network_root,
            dashboard_root,
            "cross_omics_similarity",
            "RNA-protein module similarity",
            "kme_weighted_jaccard_heatmap.png",
            "cross-omics module similarity",
            "image",
        )

        files: list[dict[str, str]] = []
        for spec in (
            ("analysis_summary", "findModules analysis summary", "05_find_modules/" + MODULE_CONDITION + "/analysis_summary.json", "summary"),
            ("module_summary", "Module summary", "05_find_modules/" + MODULE_CONDITION + "/module_summary.csv", "module"),
            ("module_assignments", "Module assignments and kME", "05_find_modules/" + MODULE_CONDITION + "/module_assignments.csv.gz", "module"),
            ("dynamic_mes", "Dynamic module eigengenes", "05_find_modules/" + MODULE_CONDITION + "/dynamic_module_eigengenes.csv", "module"),
            ("final_mes", "Final module eigengenes", "05_find_modules/" + MODULE_CONDITION + "/module_eigengenes.csv", "module"),
            ("me_linkage", "ME linkage", "05_find_modules/" + MODULE_CONDITION + "/module_eigengene_linkage.csv", "module"),
            ("soft_threshold_table", "Soft-threshold table", "05_find_modules/" + MODULE_CONDITION + "/soft_threshold_table.csv", "module"),
            ("preprocess_summary", "Selected-cut preprocess summary", "03_preprocess/02_selected_cut/preprocess_summary.json", "preprocess"),
            ("removed_samples", "Removed samples", "03_preprocess/02_selected_cut/removed_samples.csv", "preprocess"),
            ("before_outliers", "Before-preprocess outliers", "02_before_preprocess/01_descriptive_statistics/outlier_candidates.csv", "QC"),
            ("after_outliers", "After-preprocess outliers", "04_after_preprocess/01_descriptive_statistics/outlier_candidates.csv", "QC"),
            ("run_manifest", "findModules run manifest", "05_find_modules/" + MODULE_CONDITION + "/run_manifest.json", "provenance"),
            ("me_manifest", "ME dendrogram manifest", "05_find_modules/" + MODULE_CONDITION + "/me_dendrogram_manifest.json", "provenance"),
            ("official_kme", "Official all-module kME", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/feature_module_kme.csv.gz", "kME"),
            ("kme_summary", "Official kME summary", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_kme_summary.csv", "kME"),
            ("diagnosis_statistics", "One-hot diagnosis statistics (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_statistics.csv", "module diagnosis"),
            ("diagnosis_correlations", "One-hot diagnosis correlations (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_correlations.csv", "module diagnosis"),
            ("diagnosis_pvalues", "One-hot diagnosis p-values (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pvalues.csv", "module diagnosis"),
            ("diagnosis_encoding", "One-hot diagnosis encoding (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/diagnosis_trait_encoding.csv", "module diagnosis"),
            ("diagnosis_pairwise_statistics", "Pairwise disease-vs-HD statistics", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_statistics.csv", "module diagnosis"),
            ("diagnosis_pairwise_correlations", "Pairwise disease-vs-HD correlations", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_correlations.csv", "module diagnosis"),
            ("diagnosis_pairwise_pvalues", "Pairwise disease-vs-HD p-values", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_pvalues.csv", "module diagnosis"),
            ("diagnosis_pairwise_qvalues", "Pairwise disease-vs-HD q-values", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_qvalues.csv", "module diagnosis"),
            ("diagnosis_pairwise_counts", "Pairwise disease-vs-HD sample counts", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_sample_counts.csv", "module diagnosis"),
            ("diagnosis_pairwise_summary", "Pairwise disease-vs-HD analysis summary", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/pairwise_hd_analysis_summary.json", "summary"),
            ("diagnosis_pairwise_manifest", "Pairwise disease-vs-HD run manifest", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/pairwise_hd_run_manifest.json", "provenance"),
            ("diagnosis_manifest", "Official kME/one-hot diagnosis run manifest", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/run_manifest.json", "provenance"),
        ):
            existing_file(files, network_root, dashboard_root, *spec)
        for spec in (
            ("cross_pair_similarity", "All RNA-protein module pairs", "module_pair_similarity.csv", "cross-omics"),
            ("cross_best_matches", "Transcript module best matches", "transcript_module_best_matches.csv", "cross-omics"),
            ("cross_gene_contributions", "Reciprocal-best gene contributions", "reciprocal_best_gene_contributions.csv.gz", "cross-omics"),
            ("cross_mapping_audit", "Cross-omics gene mapping audit", "gene_mapping_audit.json", "cross-omics"),
            ("cross_summary", "Cross-omics analysis summary", "analysis_summary.json", "cross-omics"),
            ("cross_manifest", "Cross-omics run manifest", "run_manifest.json", "provenance"),
        ):
            existing_file(files, cross_network_root, dashboard_root, *spec)

        flags = []
        if analysis.get("collapsed_to_one_proper_module"):
            flags.append({"level": "warning", "code": "collapsed", "label": "Final modules collapsed to one"})
        unassigned_fraction = analysis.get("raw_unassigned_fraction")
        if isinstance(unassigned_fraction, (int, float)) and unassigned_fraction >= 0.5:
            flags.append({"level": "warning", "code": "high_unassigned", "label": "Raw unassigned fraction >= 50%"})
        if selected_sft and isinstance(selected_sft.get("r2"), (int, float)) and selected_sft["r2"] < 0.9:
            flags.append({"level": "note", "code": "sft_below_target", "label": "Selected power SFT R2 < 0.90"})
        if preprocess.get("finite_cut_applied"):
            flags.append({"level": "note", "code": "finite_cut", "label": "Finite sample cut applied"})
        if non_estimable_diagnosis:
            flags.append({"level": "warning", "code": "non_estimable_diagnosis", "label": f"{non_estimable_diagnosis} diagnosis tests are NA"})
        if diagnosis_q_values and min(diagnosis_q_values) < 0.05:
            flags.append({"level": "note", "code": "diagnosis_q_lt_005", "label": "At least one within-network diagnosis q < 0.05"})

        record_id = "::".join(
            [
                "transcriptomics_variance",
                row["threshold"],
                row["strategy"],
                row["mode"],
                row["specimen"],
                row["cell_type"],
            ]
        )
        records.append(
            {
                "id": record_id,
                "source_id": "transcriptomics_variance_findmodules",
                "title": f"{row['cell_type']} · {row['specimen']}",
                "dimensions": {
                    "threshold": row["threshold"],
                    "strategy": row["strategy"],
                    "mode": row["mode"],
                    "specimen": row["specimen"],
                    "cell_type": row["cell_type"],
                    "condition": analysis.get("condition"),
                },
                "status": analysis.get("status"),
                "metrics": {
                    "samples": analysis.get("samples"),
                    "genes": analysis.get("genes"),
                    "selected_power": selected_power,
                    "selected_power_sft_r2": selected_sft.get("r2") if selected_sft else None,
                    "selected_power_mean_connectivity": selected_sft.get("mean_connectivity") if selected_sft else None,
                    "dynamic_modules": analysis.get("dynamic_module_count"),
                    "final_modules": analysis.get("final_module_count"),
                    "unassigned_fraction": analysis.get("raw_unassigned_fraction"),
                    "largest_module_fraction": analysis.get("largest_module_fraction"),
                    "median_abs_kme": analysis.get("median_assigned_abs_kME"),
                    "runtime_seconds": analysis.get("runtime_seconds"),
                    "finite_cut_applied": preprocess.get("finite_cut_applied"),
                    "selected_cut": preprocess.get("selected_cut"),
                    "samples_removed_by_cut": preprocess.get("samples_removed_by_dendrogram_cut"),
                    "distance_flagged_before_cut": preprocess.get("distance_flagged_before_cut"),
                    "collapsed": analysis.get("collapsed_to_one_proper_module"),
                    "diagnosis_tests": len(diagnosis_rows) or None,
                    "minimum_diagnosis_q": min(diagnosis_q_values) if diagnosis_q_values else None,
                    "cross_omics_common_genes": cross_summary.get("common_gene_universe"),
                    "cross_omics_max_similarity": cross_summary.get("maximum_weighted_jaccard"),
                    "cross_omics_median_best_similarity": cross_summary.get("median_transcript_best_weighted_jaccard"),
                    "cross_omics_reciprocal_pairs": cross_summary.get("reciprocal_best_pairs"),
                },
                "flags": flags,
                "artifacts": artifacts,
                "files": files,
                "modules": modules,
                "provenance": {
                    "network_id": row["network_id"],
                    "analysis_finished_at_utc": analysis.get("finished_at_utc"),
                    "feature_filtering_before_findModules": analysis.get("feature_filtering_before_findModules"),
                    "expression_transform_before_findModules": analysis.get("expression_transform_before_findModules"),
                    "network_type": analysis.get("networkType"),
                    "tom_type": analysis.get("TOMType"),
                    "deep_split": analysis.get("deepSplit"),
                    "min_module_size": analysis.get("minModuleSize"),
                    "me_diss_threshold": analysis.get("MEDissThres"),
                    "cut_selection_reason": preprocess.get("cut_selection_reason"),
                    "kme_method": diagnosis_summary.get("kme_method"),
                    "module_diagnosis_method": pairwise_summary.get("method") or diagnosis_summary.get("heatmap_method"),
                    "cross_omics_method": "absolute own-module kME weighted Jaccard" if cross_summary else None,
                    "cross_omics_proteomics_condition": cross_summary.get("proteomics_condition"),
                },
            }
        )

    source = {
        "id": "transcriptomics_variance_findmodules",
        "label": "Transcriptomics variance-filter findModules",
        "stage": "04_variance_filtered",
        "summary_path": dashboard_relative(summary_path, dashboard_root),
        "record_count": len(records),
    }
    return source, records


def build_proteomics_modules(
    run_root: Path,
    dashboard_root: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    stage_root = run_root / "experiments/01_batch_corrected/proteomics"
    records: list[dict[str, Any]] = []
    for specimen in ("CSF", "PB"):
        condition_rows = {
            row["condition"]: row
            for row in read_csv_rows(stage_root / "05_find_modules" / specimen / "condition_summary.csv")
        }
        for condition, condition_row in condition_rows.items():
            module_root = stage_root / "05_find_modules" / specimen / condition
            diagnosis_root = stage_root / "06_kme_module_diagnosis" / specimen / condition
            analysis = read_json(module_root / "analysis_summary.json")
            diagnosis_summary = read_json(diagnosis_root / "analysis_summary.json")
            pairwise_summary = read_json(diagnosis_root / "pairwise_hd_analysis_summary.json")
            preprocess_path = stage_root / "03_preprocess" / specimen / "02_selected_cut" / "preprocess_summary.json"
            preprocess = read_json(preprocess_path) if preprocess_path.exists() else {}
            diagnosis_rows = read_csv_rows(
                diagnosis_root / "module_diagnosis_pairwise_hd_statistics.csv"
            )
            diagnosis_q_values = [
                value
                for value in (
                    parse_number(item.get("bh_fdr_q_within_network")) for item in diagnosis_rows
                )
                if value is not None
            ]
            non_estimable_diagnosis = sum(
                parse_number(item.get("p_value_two_sided")) is None for item in diagnosis_rows
            )

            modules = []
            for module in read_csv_rows(module_root / "module_summary.csv"):
                modules.append(
                    {
                        "color": module.get("module_color"),
                        "features": parse_number(module.get("feature_count")),
                        "assigned_features": parse_number(module.get("assigned_feature_count")),
                        "unassigned_features": parse_number(module.get("raw_unassigned_feature_count")),
                        "median_abs_kme": parse_number(
                            module.get("assigned_abs_kME_median")
                            or module.get("assigned_kme_median")
                        ),
                        "fully_unassigned": str(module.get("is_fully_raw_unassigned", "")).lower() == "true",
                    }
                )

            selected_power = parse_number(analysis.get("selected_power"))
            selected_sft = None
            for metric in read_csv_rows(module_root / "soft_threshold_table.csv"):
                if parse_number(metric.get("Power")) == selected_power:
                    selected_sft = {
                        "r2": parse_number(metric.get("SFT.R.sq")),
                        "mean_connectivity": parse_number(metric.get("mean(k)")),
                    }
                    break

            artifacts: list[dict[str, str]] = []
            mirrored_image_artifact(
                artifacts,
                diagnosis_root / "figures/module_diagnosis_pairwise_hd_heatmap.png",
                dashboard_root,
                "diagnosis_pairwise_hd",
                "Module-diagnosis (vs HD)",
                "module diagnosis",
            )
            for spec in (
                ("diagnosis_heatmap", "One-hot diagnosis (audit)", f"06_kme_module_diagnosis/{specimen}/{condition}/figures/module_diagnosis_heatmap.png", "module diagnosis", "image"),
                ("me_dendrogram", "ME dendrogram", f"05_find_modules/{specimen}/{condition}/figures/module_eigengene_dendrogram.png", "findModules", "image"),
                ("gene_dendrogram", "Gene dendrogram", f"05_find_modules/{specimen}/{condition}/figures/gene_dendrogram_modules.png", "findModules", "image"),
                ("soft_threshold", "Soft threshold", f"05_find_modules/{specimen}/{condition}/figures/soft_threshold.png", "findModules", "image"),
                ("module_sizes", "Module sizes", f"05_find_modules/{specimen}/{condition}/figures/module_sizes.png", "findModules", "image"),
                ("selected_cut", "Selected cut", f"03_preprocess/{specimen}/02_selected_cut/dendrogram_selected_cut.png", "preprocess", "image"),
                ("after_heatmap", "After heatmap", f"04_after_preprocess/{specimen}/02_data_distribution/heatmap.png", "after preprocess", "image"),
                ("after_pca", "After PCA", f"04_after_preprocess/{specimen}/02_data_distribution/pca_metadata.png", "after preprocess", "image"),
                ("after_pca_3d", "After PCA 3D", f"04_after_preprocess/{specimen}/02_data_distribution/pca_3d_interactive.html", "after preprocess", "html"),
            ):
                existing_artifact(artifacts, stage_root, dashboard_root, *spec)

            files: list[dict[str, str]] = []
            for spec in (
                ("analysis_summary", "findModules analysis summary", f"05_find_modules/{specimen}/{condition}/analysis_summary.json", "summary"),
                ("module_summary", "Module summary", f"05_find_modules/{specimen}/{condition}/module_summary.csv", "module"),
                ("source_kme", "findModules feature-module kME", f"05_find_modules/{specimen}/{condition}/feature_module_kme.csv", "module"),
                ("final_mes", "Final module eigengenes", f"05_find_modules/{specimen}/{condition}/module_eigengenes.csv", "module"),
                ("official_kme", "Official rerun all-module kME", f"06_kme_module_diagnosis/{specimen}/{condition}/feature_module_kme.csv.gz", "kME"),
                ("kme_summary", "Official kME summary", f"06_kme_module_diagnosis/{specimen}/{condition}/module_kme_summary.csv", "kME"),
                ("diagnosis_statistics", "One-hot diagnosis statistics (audit)", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_statistics.csv", "module diagnosis"),
                ("diagnosis_correlations", "One-hot diagnosis correlations (audit)", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_correlations.csv", "module diagnosis"),
                ("diagnosis_pvalues", "One-hot diagnosis p-values (audit)", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pvalues.csv", "module diagnosis"),
                ("diagnosis_encoding", "One-hot diagnosis encoding (audit)", f"06_kme_module_diagnosis/{specimen}/{condition}/diagnosis_trait_encoding.csv", "module diagnosis"),
                ("diagnosis_pairwise_statistics", "Pairwise disease-vs-HD statistics", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pairwise_hd_statistics.csv", "module diagnosis"),
                ("diagnosis_pairwise_correlations", "Pairwise disease-vs-HD correlations", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pairwise_hd_correlations.csv", "module diagnosis"),
                ("diagnosis_pairwise_pvalues", "Pairwise disease-vs-HD p-values", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pairwise_hd_pvalues.csv", "module diagnosis"),
                ("diagnosis_pairwise_qvalues", "Pairwise disease-vs-HD q-values", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pairwise_hd_qvalues.csv", "module diagnosis"),
                ("diagnosis_pairwise_counts", "Pairwise disease-vs-HD sample counts", f"06_kme_module_diagnosis/{specimen}/{condition}/module_diagnosis_pairwise_hd_sample_counts.csv", "module diagnosis"),
                ("diagnosis_pairwise_summary", "Pairwise disease-vs-HD analysis summary", f"06_kme_module_diagnosis/{specimen}/{condition}/pairwise_hd_analysis_summary.json", "summary"),
                ("diagnosis_pairwise_manifest", "Pairwise disease-vs-HD run manifest", f"06_kme_module_diagnosis/{specimen}/{condition}/pairwise_hd_run_manifest.json", "provenance"),
                ("findmodules_manifest", "findModules run manifest", f"05_find_modules/{specimen}/{condition}/run_manifest.json", "provenance"),
                ("diagnosis_manifest", "Official kME/one-hot diagnosis run manifest", f"06_kme_module_diagnosis/{specimen}/{condition}/run_manifest.json", "provenance"),
            ):
                existing_file(files, stage_root, dashboard_root, *spec)

            flags = []
            if str(condition_row.get("is_pywgcna_default_condition", "")).lower() == "true":
                flags.append({"level": "note", "code": "default_condition", "label": "PyWGCNA default module condition"})
            if condition.startswith("04_"):
                flags.append({"level": "note", "code": "reproduction_target", "label": "Colleague-result reproduction target"})
            if analysis.get("collapsed_to_one_proper_module"):
                flags.append({"level": "warning", "code": "collapsed", "label": "Final modules collapsed to one"})
            if non_estimable_diagnosis:
                flags.append({"level": "warning", "code": "non_estimable_diagnosis", "label": f"{non_estimable_diagnosis} diagnosis tests are NA"})
            if diagnosis_q_values and min(diagnosis_q_values) < 0.05:
                flags.append({"level": "note", "code": "diagnosis_q_lt_005", "label": "At least one within-network diagnosis q < 0.05"})

            records.append(
                {
                    "id": f"proteomics::{specimen}::{condition}",
                    "source_id": "proteomics_module_conditions",
                    "title": f"Proteomics {specimen} · {condition}",
                    "dimensions": {
                        "threshold": "",
                        "strategy": condition,
                        "mode": "proteomics",
                        "specimen": specimen,
                        "cell_type": "ALL",
                        "condition": condition,
                    },
                    "status": analysis.get("status"),
                    "metrics": {
                        "samples": analysis.get("samples"),
                        "genes": analysis.get("features_used"),
                        "selected_power": selected_power,
                        "selected_power_sft_r2": selected_sft.get("r2") if selected_sft else None,
                        "selected_power_mean_connectivity": selected_sft.get("mean_connectivity") if selected_sft else None,
                        "dynamic_modules": analysis.get("dynamic_module_count"),
                        "final_modules": analysis.get("final_module_count"),
                        "unassigned_fraction": analysis.get("raw_unassigned_fraction"),
                        "largest_module_fraction": analysis.get("largest_module_fraction"),
                        "median_abs_kme": analysis.get("median_assigned_abs_kME"),
                        "runtime_seconds": None,
                        "finite_cut_applied": True,
                        "selected_cut": preprocess.get("selected_cut"),
                        "samples_removed_by_cut": preprocess.get("samples_removed_by_dendrogram_cut"),
                        "distance_flagged_before_cut": preprocess.get("distance_flagged_before_cut"),
                        "collapsed": analysis.get("collapsed_to_one_proper_module"),
                        "diagnosis_tests": len(diagnosis_rows),
                        "minimum_diagnosis_q": min(diagnosis_q_values) if diagnosis_q_values else None,
                    },
                    "flags": flags,
                    "artifacts": artifacts,
                    "files": files,
                    "modules": modules,
                    "provenance": {
                        "network_id": f"proteomics_{specimen}_{condition}",
                        "analysis_finished_at_utc": analysis.get("finished_at_utc"),
                        "feature_filtering_before_findModules": analysis.get("feature_filtering_before_findModules"),
                        "network_type": analysis.get("networkType"),
                        "tom_type": analysis.get("TOMType"),
                        "deep_split": analysis.get("deepSplit"),
                        "min_module_size": analysis.get("minModuleSize"),
                        "me_diss_threshold": analysis.get("MEDissThres"),
                        "kme_method": diagnosis_summary.get("kme_method"),
                        "module_diagnosis_method": pairwise_summary.get("method") or diagnosis_summary.get("heatmap_method"),
                    },
                }
            )

    source = {
        "id": "proteomics_module_conditions",
        "label": "Proteomics module parameter conditions",
        "stage": "01_batch_corrected",
        "summary_path": dashboard_relative(
            stage_root / "06_kme_module_diagnosis" / "summary.csv", dashboard_root
        ),
        "record_count": len(records),
    }
    return source, records


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_records(records: list[dict[str, Any]], dashboard_root: Path) -> None:
    """Fail before publishing a dashboard with ambiguous or broken records."""

    record_ids = [record["id"] for record in records]
    duplicated = sorted({item for item in record_ids if record_ids.count(item) > 1})
    if duplicated:
        raise ValueError(f"Duplicate dashboard record IDs: {duplicated[:5]}")

    missing_paths: list[str] = []
    missing_me: list[str] = []
    for record in records:
        items = [*record.get("artifacts", []), *record.get("files", [])]
        for item in items:
            if not (dashboard_root / item["path"]).resolve().exists():
                missing_paths.append(f"{record['id']}::{item['id']}::{item['path']}")
        if record.get("status") == "complete_findModules" and not any(
            item["id"] == "me_dendrogram" for item in record.get("artifacts", [])
        ):
            missing_me.append(record["id"])

    if missing_paths:
        raise FileNotFoundError(
            "Dashboard contains missing artifact/file links: "
            + "; ".join(missing_paths[:10])
        )
    if missing_me:
        raise ValueError(
            "Completed findModules records without an ME dendrogram: "
            + ", ".join(missing_me[:10])
        )


def build(run_root: Path) -> None:
    dashboard_root = run_root / "dashboard"
    dashboard_root.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).isoformat()

    # Add future experiment families here as another adapter returning the same schema.
    transcriptomics_source, transcriptomics_records = build_variance_filtered_transcriptomics(
        run_root, dashboard_root
    )
    proteomics_source, proteomics_records = build_proteomics_modules(run_root, dashboard_root)
    sources = [transcriptomics_source, proteomics_source]
    records = [*transcriptomics_records, *proteomics_records]
    validate_records(records, dashboard_root)
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "run_root": run_root.as_posix(),
        "sources": sources,
        "records": records,
    }
    data_path = dashboard_root / "dashboard-data.js"
    data_path.write_text(
        "window.EXPERIMENT_DASHBOARD_DATA = "
        + json.dumps(bundle, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )
    html_path = dashboard_root / "index.html"
    dashboard_html = DASHBOARD_HTML.replace(
        "(artifacts.find(a=>a.id==='diagnosis_heatmap')||artifacts.find(a=>a.id==='me_dendrogram')||artifacts[0])?.id",
        "(artifacts.find(a=>a.id==='diagnosis_pairwise_hd')||artifacts.find(a=>a.id==='cross_omics_similarity')||artifacts.find(a=>a.id==='diagnosis_heatmap')||artifacts.find(a=>a.id==='me_dendrogram')||artifacts[0])?.id",
    ).replace(
        "['Samples cut',fmt(r.metrics.samples_removed_by_cut)]",
        "['Common omics genes',fmt(r.metrics.cross_omics_common_genes)],['Max weighted Jaccard',fmt(r.metrics.cross_omics_max_similarity,3)],['Median best match',fmt(r.metrics.cross_omics_median_best_similarity,3)],['Reciprocal pairs',fmt(r.metrics.cross_omics_reciprocal_pairs)],['Samples cut',fmt(r.metrics.samples_removed_by_cut)]",
    ).replace(
        "</style>",
        DESIGN_SYSTEM_CSS + "\n</style>",
        1,
    )
    html_path.write_text(dashboard_html, encoding="utf-8")
    readme_path = dashboard_root / "README.md"
    readme_path.write_text(
        """# Experiment inspection dashboard

`index.html`을 직접 열면 된다. 별도 web server나 CDN은 필요하지 않는다.

## 구성

- `index.html`: 바뀌지 않는 inspection UI
- `dashboard-data.js`: 각 실험 결과를 공통 schema로 정규화한 재생성 가능 데이터
- `dashboard_manifest.json`: 생성 시각, record 수, source와 builder hash
- `dashboard_preview.png`: 브라우저 렌더링 확인용 예시 화면

현재 예시는 transcriptomics variance-filter `findModules` 168개 record와
proteomics module-parameter condition 8개 record를 연결한다.
기본 **Cell type comparison** 화면에서는 cell type을 먼저 고르고, 산출물 탭을
바꾸면 모든 threshold/feature-selection/specimen 조건의 같은 산출물이 한 화면에
정렬된다. **Individual run inspection**으로 전환하면 기존의 단일 record 검토,
module table과 provenance를 그대로 사용할 수 있다.
판단, 체크리스트와 메모는 브라우저 `localStorage`에 자동 저장되며
**Export review CSV**로 이식 가능한 의사결정 기록을 남긴다. 다른 브라우저나
컴퓨터에서는 그 CSV를 **Import review CSV**로 다시 불러올 수 있다.

## 갱신

결과 파일이 추가되거나 바뀌면 다음 명령으로 데이터 bundle과 manifest를 갱신한다.

```powershell
uv run python results_20261119/code/src/build_experiment_inspection_dashboard.py --run-root results_20261119
```

새 실험 family를 포함하려면 builder에 source adapter를 하나 추가해 동일한
record/artifact/file/provenance schema를 반환하게 한다. UI는 수정하지 않는다.

의사결정의 원본은 export한 CSV다. `localStorage`만 사용한 판단은 해당 브라우저에만
남으므로 검토가 끝난 시점에는 반드시 CSV를 export한다.
""",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "records": len(records),
        "sources": sources,
        "outputs": [html_path.name, data_path.name, readme_path.name],
        "builder": {
            "path": Path(__file__).resolve().as_posix(),
            "sha256": sha256(Path(__file__).resolve()),
        },
    }
    (dashboard_root / "dashboard_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Built {html_path} with {len(records)} experiment records")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    return parser.parse_args()


DASHBOARD_HTML = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Experiment inspection dashboard</title>
<style>
:root{--ink:#18212f;--muted:#667085;--line:#d9dee7;--panel:#f5f7fa;--accent:#0f6f78;--accent2:#dceff0;--warn:#a54b12;--danger:#9b2c2c;--ok:#276749;--white:#fff}
*{box-sizing:border-box}body{margin:0;font:14px/1.45 Arial,Helvetica,sans-serif;color:var(--ink);background:var(--white)}button,select,input,textarea{font:inherit}button,select,input,textarea{border:1px solid var(--line);background:white;color:var(--ink);border-radius:3px}button{cursor:pointer}.top{padding:18px 24px 14px;border-bottom:1px solid var(--line)}h1{font-size:24px;margin:0 0 4px}.sub{color:var(--muted)}.kpis{display:grid;grid-template-columns:repeat(6,minmax(110px,1fr));gap:1px;background:var(--line);border-bottom:1px solid var(--line)}.kpi{background:white;padding:12px 18px}.kpi b{display:block;font-size:21px}.kpi span{color:var(--muted);font-size:12px}.toolbar{position:sticky;top:0;z-index:10;display:flex;gap:8px;align-items:center;flex-wrap:wrap;padding:10px 16px;background:var(--panel);border-bottom:1px solid var(--line)}.toolbar select,.toolbar input{padding:7px 8px;min-width:130px}.toolbar #search{min-width:230px}.toolbar button{padding:7px 10px}.spacer{flex:1}.layout{display:grid;grid-template-columns:360px minmax(0,1fr);min-height:calc(100vh - 210px)}.list{border-right:1px solid var(--line);height:calc(100vh - 210px);overflow:auto;background:#fbfcfd}.list-head{padding:10px 14px;color:var(--muted);border-bottom:1px solid var(--line);position:sticky;top:0;background:#fbfcfd}.record{width:100%;text-align:left;padding:12px 14px;border:0;border-bottom:1px solid var(--line);background:transparent;border-radius:0}.record:hover,.record.active{background:var(--accent2)}.record-title{font-weight:700}.record-meta,.record-metrics{display:flex;gap:8px;flex-wrap:wrap;color:var(--muted);font-size:12px;margin-top:4px}.pill{display:inline-block;padding:1px 6px;border:1px solid var(--line);border-radius:12px}.pill.warn{color:var(--warn);border-color:#e7b58f}.pill.ok{color:var(--ok)}.detail{padding:20px 24px 50px;min-width:0}.detail-head{display:flex;gap:20px;align-items:flex-start}.detail-head h2{margin:0 0 6px;font-size:22px}.detail-head .actions{margin-left:auto;display:flex;gap:8px}.detail-head button{padding:7px 10px}.metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(135px,1fr));gap:1px;background:var(--line);border:1px solid var(--line);margin:16px 0}.metric{background:white;padding:10px 12px}.metric .label{font-size:11px;color:var(--muted)}.metric .value{font-size:16px;font-weight:700}.flags{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0 14px}.flag{padding:4px 8px;background:#fff7ed;color:var(--warn);border:1px solid #fed7aa}.review{display:grid;grid-template-columns:minmax(310px,1fr) minmax(320px,1fr);gap:18px;padding:14px;background:var(--panel);border:1px solid var(--line);margin:14px 0}.decision-row{display:flex;gap:7px;margin:8px 0}.decision-row button{padding:7px 10px}.decision-row button.selected{background:var(--accent);color:white;border-color:var(--accent)}.checklist label{display:block;margin:5px 0}.review textarea{width:100%;min-height:88px;padding:8px;resize:vertical}.section-title{font-size:16px;margin:22px 0 8px}.artifact-tabs{display:flex;gap:6px;overflow:auto;padding-bottom:6px}.artifact-tabs button{padding:6px 9px;white-space:nowrap}.artifact-tabs button.active{background:var(--accent);color:white;border-color:var(--accent)}.viewer{border:1px solid var(--line);min-height:420px;background:white;display:flex;align-items:flex-start;justify-content:center;overflow:auto}.viewer img{display:block;max-width:100%;height:auto}.viewer iframe{width:100%;height:680px;border:0}.viewer .open{padding:18px}.lower{display:grid;grid-template-columns:minmax(0,1.2fr) minmax(300px,.8fr);gap:18px}.table-wrap{overflow:auto;border:1px solid var(--line);max-height:360px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}th{position:sticky;top:0;background:var(--panel)}.swatch{display:inline-block;width:11px;height:11px;border:1px solid #777;margin-right:6px;vertical-align:-1px}.files{border:1px solid var(--line);padding:10px 12px}.files details{border-bottom:1px solid var(--line);padding:7px 0}.files details:last-child{border-bottom:0}.files a{display:block;color:var(--accent);padding:3px 0;word-break:break-all}.compare{display:none;margin:0 16px 20px;border:1px solid var(--line);padding:14px}.compare.visible{display:block}.compare-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}.compare-card{border:1px solid var(--line);padding:10px}.compare-card img{width:100%;height:auto}.compare-card h3{font-size:14px;margin:0 0 7px}.empty{padding:50px;color:var(--muted);text-align:center}@media(max-width:1000px){.kpis{grid-template-columns:repeat(3,1fr)}.layout{grid-template-columns:1fr}.list{height:320px;border-right:0;border-bottom:1px solid var(--line)}.review,.lower{grid-template-columns:1fr}}
</style>
</head>
<body data-view="cell">
<header class="top"><h1>Experiment inspection dashboard</h1><div class="sub" id="sourceLine"></div></header>
<section class="kpis" id="kpis"></section>
<div class="toolbar">
  <select id="viewMode"><option value="cell">Cell type comparison</option><option value="run">Individual run inspection</option></select><select id="source"><option value="">All sources</option></select><select id="cell"><option value="">Choose cell type</option></select><select class="run-only" id="threshold"><option value="">All thresholds</option></select><select class="run-only" id="strategy"><option value="">All strategies</option></select><select class="run-only" id="mode"><option value="">All modes</option></select><select class="run-only" id="specimen"><option value="">All specimens</option></select><select class="run-only" id="status"><option value="">All statuses</option></select><select class="run-only" id="decision"><option value="">All decisions</option><option value="unreviewed">Unreviewed</option><option value="use">Use</option><option value="hold">Hold</option><option value="reject">Do not use</option></select><select class="run-only" id="sort"><option value="default">Default order</option><option value="modules_desc">Final modules ↓</option><option value="unassigned_desc">Unassigned ↓</option><option value="runtime_desc">Runtime ↓</option></select><input class="run-only" id="search" placeholder="Search cell type or network"><span class="spacer"></span><button class="run-only" id="nextUnreviewed">Next unreviewed</button><button class="run-only" id="showCompare">Compare <span id="compareCount">0</span></button><button id="exportReview">Export review CSV</button><label><button type="button" id="importButton">Import review CSV</button><input id="importReview" type="file" accept=".csv" hidden></label>
</div>
<main class="cell-comparison cell-only" id="cellComparison"></main>
<section class="compare run-only" id="comparePanel"><div class="detail-head"><h2>Comparison workspace</h2><div class="actions"><button id="clearCompare">Clear</button></div></div><div class="compare-grid" id="compareGrid"></div></section>
<main class="layout run-only" id="runInspection"><aside class="list"><div class="list-head" id="listCount"></div><div id="recordList"></div></aside><section class="detail" id="detail"></section></main>
<script src="dashboard-data.js"></script>
<script>
const bundle=window.EXPERIMENT_DASHBOARD_DATA;const records=bundle.records;const storageKey=`experiment-review:${bundle.schema_version}:${bundle.run_root}`;let reviews=JSON.parse(localStorage.getItem(storageKey)||'{}');let selectedId=records[0]?.id||null;let artifactChoice={};let comparisonArtifactChoice=null;let compareIds=[];
const $=id=>document.getElementById(id);const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const fmt=(v,d=2)=>v===null||v===undefined||v===''?'—':typeof v==='number'?v.toLocaleString(undefined,{maximumFractionDigits:d}):esc(v);const stat=v=>v===null||v===undefined||v===''?'—':typeof v==='number'&&Math.abs(v)>0&&Math.abs(v)<0.0001?v.toExponential(2):fmt(v,4);const pct=v=>v===null||v===undefined?'—':`${(v*100).toFixed(1)}%`;
const artifactOrder=['feature_selection','before_heatmap','before_pca','before_pca_3d','cut_diagnostic','selected_cut','after_heatmap','after_pca','after_pca_3d','soft_threshold','gene_dendrogram','me_dendrogram','module_sizes','diagnosis_pairwise_hd','diagnosis_heatmap','cross_omics_similarity'];
function reviewOf(id){return reviews[id]||{decision:'unreviewed',notes:'',checks:{power:false,modules:false,outliers:false,ready:false},updated_at:null}}
function saveReview(id,patch){reviews[id]={...reviewOf(id),...patch,updated_at:new Date().toISOString()};localStorage.setItem(storageKey,JSON.stringify(reviews));renderAll()}
function unique(path){return [...new Set(records.map(r=>path(r)).filter(Boolean))].sort()}
function fillSelect(id,values){const node=$(id);for(const v of values){const o=document.createElement('option');o.value=v;o.textContent=v;node.appendChild(o)}node.addEventListener('change',renderAll)}
fillSelect('source',bundle.sources.map(s=>s.id));fillSelect('threshold',unique(r=>r.dimensions.threshold));fillSelect('strategy',unique(r=>r.dimensions.strategy));fillSelect('mode',unique(r=>r.dimensions.mode));fillSelect('specimen',unique(r=>r.dimensions.specimen));fillSelect('cell',unique(r=>r.dimensions.cell_type));fillSelect('status',unique(r=>r.status));for(const option of $('source').options){const source=bundle.sources.find(s=>s.id===option.value);if(source)option.textContent=source.label}for(const option of $('cell').options)if(option.value)option.textContent=prettyDimension(option.value);['decision','sort'].forEach(id=>$(id).addEventListener('change',renderAll));$('search').addEventListener('input',renderAll);
$('viewMode').addEventListener('change',renderAll);
const transcriptomicsSource=bundle.sources.find(s=>s.id.includes('transcriptomics'))?.id||bundle.sources[0]?.id||'';$('source').value=transcriptomicsSource;if([...$('cell').options].some(o=>o.value==='B_cell'))$('cell').value='B_cell';
function ensureCellForSource(){const source=$('source').value,available=unique(r=>(!source||r.source_id===source)?r.dimensions.cell_type:null);if(!available.includes($('cell').value))$('cell').value=available.includes('B_cell')?'B_cell':available[0]||''}
$('source').addEventListener('change',()=>{if($('viewMode').value==='cell'){ensureCellForSource();renderAll()}});
function filtered(){const q=$('search').value.trim().toLowerCase();let out=records.filter(r=>(!$('source').value||r.source_id===$('source').value)&&(!$('threshold').value||r.dimensions.threshold===$('threshold').value)&&(!$('strategy').value||r.dimensions.strategy===$('strategy').value)&&(!$('mode').value||r.dimensions.mode===$('mode').value)&&(!$('specimen').value||r.dimensions.specimen===$('specimen').value)&&(!$('cell').value||r.dimensions.cell_type===$('cell').value)&&(!$('status').value||r.status===$('status').value)&&(!$('decision').value||reviewOf(r.id).decision===$('decision').value)&&(!q||`${r.title} ${r.provenance.network_id} ${r.dimensions.strategy}`.toLowerCase().includes(q)));const s=$('sort').value;if(s==='modules_desc')out.sort((a,b)=>(b.metrics.final_modules??-1)-(a.metrics.final_modules??-1));if(s==='unassigned_desc')out.sort((a,b)=>(b.metrics.unassigned_fraction??-1)-(a.metrics.unassigned_fraction??-1));if(s==='runtime_desc')out.sort((a,b)=>(b.metrics.runtime_seconds??-1)-(a.metrics.runtime_seconds??-1));return out}
function prettyDimension(v){return String(v??'—').replace(/^\d+_/,'').replaceAll('_',' ')}
function comparisonRows(){
  const source=$('source').value,cell=$('cell').value;
  return records.filter(r=>(!source||r.source_id===source)&&cell&&r.dimensions.cell_type===cell).sort((a,b)=>{
    const numberIn=v=>Number(String(v??'').match(/\d+/)?.[0]??999),specimenRank=v=>v==='CSF'?0:v==='PB'?1:2;
    const ka=[numberIn(a.dimensions.threshold),numberIn(a.dimensions.strategy),String(a.dimensions.mode??''),specimenRank(a.dimensions.specimen),String(a.provenance.condition??'')];
    const kb=[numberIn(b.dimensions.threshold),numberIn(b.dimensions.strategy),String(b.dimensions.mode??''),specimenRank(b.dimensions.specimen),String(b.provenance.condition??'')];
    for(let i=0;i<ka.length;i++){const d=typeof ka[i]==='number'?ka[i]-kb[i]:ka[i].localeCompare(kb[i]);if(d)return d}return 0;
  });
}
function comparisonArtifacts(rows){
  const byId=new Map();for(const r of rows)for(const a of r.artifacts)if(!byId.has(a.id))byId.set(a.id,a);
  return [...byId.values()].sort((a,b)=>{const ai=artifactOrder.indexOf(a.id),bi=artifactOrder.indexOf(b.id);return (ai<0?999:ai)-(bi<0?999:bi)||a.label.localeCompare(b.label)});
}
function renderCellComparison(){
  const rows=comparisonRows(),artifacts=comparisonArtifacts(rows),cell=$('cell').value;
  if(!cell){$('cellComparison').innerHTML='<div class="empty">Choose one cell type to compare its conditions.</div>';return}
  if(!rows.length){$('cellComparison').innerHTML='<div class="empty">No experiment conditions are available for this source and cell type.</div>';return}
  if(!comparisonArtifactChoice||!artifacts.some(a=>a.id===comparisonArtifactChoice))comparisonArtifactChoice=artifacts[0]?.id||null;
  const chosenMeta=artifacts.find(a=>a.id===comparisonArtifactChoice),available=rows.filter(r=>r.artifacts.some(a=>a.id===comparisonArtifactChoice)).length;
  const tabs=artifacts.map(a=>`<button data-comparison-artifact="${esc(a.id)}" class="${a.id===comparisonArtifactChoice?'active':''}">${esc(a.label)}</button>`).join('');
  const panels=rows.map(r=>{
    const a=r.artifacts.find(x=>x.id===comparisonArtifactChoice),rv=reviewOf(r.id);
    const condition=[r.dimensions.threshold,r.dimensions.strategy,r.dimensions.specimen].filter(Boolean).map(prettyDimension).join(' · ');
    const secondary=[r.dimensions.mode,`n=${fmt(r.metrics.samples)}`,`genes=${fmt(r.metrics.genes)}`,`modules=${fmt(r.metrics.dynamic_modules)}→${fmt(r.metrics.final_modules)}`].filter(Boolean).map(prettyDimension).join(' · ');
    const visual=!a?'<div class="condition-missing">This artifact was not generated for this condition.</div>':a.kind==='html'?`<iframe loading="lazy" src="${esc(a.path)}" title="${esc(a.label)} · ${esc(condition)}"></iframe>`:`<a href="${esc(a.path)}" target="_blank"><img loading="lazy" src="${esc(a.path)}" alt="${esc(a.label)} · ${esc(condition)}"></a>`;
    const open=a?`<a href="${esc(a.path)}" target="_blank">Open artifact</a>`:'';
    return `<article class="condition-panel"><div class="condition-head"><div><h3>${esc(condition||r.title)}</h3><div class="condition-meta">${secondary}</div></div><span class="condition-review ${esc(rv.decision)}">${esc(rv.decision)}</span></div><div class="condition-viewer">${visual}</div><div class="condition-links">${open}<button data-open-record="${esc(r.id)}">Inspect this run</button></div></article>`;
  }).join('');
  $('cellComparison').innerHTML=`<div class="comparison-head"><div><h2>${esc(prettyDimension(cell))}</h2><div class="sub">The same output is aligned across every saved condition for this cell type.</div></div><div class="comparison-counts">${rows.length} conditions · ${artifacts.length} output tabs · ${available}/${rows.length} shown in this tab</div></div><div class="artifact-tabs">${tabs}</div><div class="comparison-note">Current output: <b>${esc(chosenMeta?.label||'—')}</b> · sorted by threshold → feature strategy → specimen, so matched CSF/PB conditions stay adjacent.</div><section class="condition-grid">${panels}</section>`;
  document.querySelectorAll('[data-comparison-artifact]').forEach(n=>n.onclick=()=>{comparisonArtifactChoice=n.dataset.comparisonArtifact;renderCellComparison()});
  document.querySelectorAll('[data-open-record]').forEach(n=>n.onclick=()=>{selectedId=n.dataset.openRecord;$('viewMode').value='run';renderAll()});
}
function renderKpis(){const complete=records.filter(r=>r.status==='complete_findModules').length,collapsed=records.filter(r=>r.metrics.collapsed).length,reviewed=records.filter(r=>reviewOf(r.id).decision!=='unreviewed').length,holds=records.filter(r=>reviewOf(r.id).decision==='hold').length,uses=records.filter(r=>reviewOf(r.id).decision==='use').length;const values=[[records.length,'Experiment records'],[complete,'findModules complete'],[collapsed,'Collapsed'],[reviewed,'Reviewed'],[uses,'Use'],[holds,'Hold']];$('kpis').innerHTML=values.map(([v,l])=>`<div class="kpi"><b>${v}</b><span>${l}</span></div>`).join('')}
function renderList(){const rows=filtered();$('listCount').textContent=`${rows.length} / ${records.length} records`;if(rows.length&&!rows.some(r=>r.id===selectedId))selectedId=rows[0].id;$('recordList').innerHTML=rows.map(r=>{const rv=reviewOf(r.id),warn=r.flags.some(f=>f.level==='warning');return `<button class="record ${r.id===selectedId?'active':''}" data-id="${esc(r.id)}"><div class="record-title">${esc(r.title)} ${warn?'<span class="pill warn">flag</span>':''}</div><div class="record-meta"><span>${esc(r.dimensions.threshold)}</span><span>${esc(r.dimensions.strategy)}</span><span>${esc(r.dimensions.mode)}</span></div><div class="record-metrics"><span>n=${fmt(r.metrics.samples)}</span><span>genes=${fmt(r.metrics.genes)}</span><span>dynamic=${fmt(r.metrics.dynamic_modules)}</span><span>final=${fmt(r.metrics.final_modules)}</span><span class="pill ${rv.decision==='use'?'ok':rv.decision==='hold'||rv.decision==='reject'?'warn':''}">${esc(rv.decision)}</span></div></button>`}).join('')||'<div class="empty">No matching experiments</div>';document.querySelectorAll('.record').forEach(n=>n.onclick=()=>{selectedId=n.dataset.id;renderAll()})}
function metric(label,value){return `<div class="metric"><div class="label">${esc(label)}</div><div class="value">${value}</div></div>`}
function renderDetail(){const r=records.find(x=>x.id===selectedId);if(!r){$('detail').innerHTML='<div class="empty">Select an experiment</div>';return}const rv=reviewOf(r.id),artifacts=r.artifacts;if(!artifactChoice[r.id]||!artifacts.some(a=>a.id===artifactChoice[r.id]))artifactChoice[r.id]=(artifacts.find(a=>a.id==='diagnosis_heatmap')||artifacts.find(a=>a.id==='me_dendrogram')||artifacts[0])?.id;const chosen=artifacts.find(a=>a.id===artifactChoice[r.id]);const metrics=[['Samples',fmt(r.metrics.samples)],['Genes',fmt(r.metrics.genes)],['Power',fmt(r.metrics.selected_power)],['Selected SFT R²',fmt(r.metrics.selected_power_sft_r2,3)],['Mean connectivity',fmt(r.metrics.selected_power_mean_connectivity,1)],['Dynamic modules',fmt(r.metrics.dynamic_modules)],['Final modules',fmt(r.metrics.final_modules)],['Unassigned',pct(r.metrics.unassigned_fraction)],['Largest module',pct(r.metrics.largest_module_fraction)],['Median |kME|',fmt(r.metrics.median_abs_kme,3)],['Diagnosis tests',fmt(r.metrics.diagnosis_tests)],['Min diagnosis q',stat(r.metrics.minimum_diagnosis_q)],['Samples cut',fmt(r.metrics.samples_removed_by_cut)],['Runtime (s)',fmt(r.metrics.runtime_seconds,1)]];const groups={};for(const f of r.files)(groups[f.category]??=[]).push(f);const proper=r.modules.filter(m=>!m.fully_unassigned).sort((a,b)=>(b.features??0)-(a.features??0));$('detail').innerHTML=`<div class="detail-head"><div><h2>${esc(r.title)}</h2><div class="sub">${esc(r.provenance.network_id)} · ${esc(r.status)}</div></div><div class="actions"><button id="addCompare">${compareIds.includes(r.id)?'Remove from compare':'Add to compare'}</button></div></div><div class="flags">${r.flags.map(f=>`<span class="flag">${esc(f.label)}</span>`).join('')}</div><div class="metric-grid">${metrics.map(x=>metric(...x)).join('')}</div><section class="review"><div><b>Decision</b><div class="decision-row">${[['use','Use'],['hold','Hold'],['reject','Do not use'],['unreviewed','Clear']].map(([v,l])=>`<button data-decision="${v}" class="${rv.decision===v?'selected':''}">${l}</button>`).join('')}</div><div class="checklist"><label><input type="checkbox" data-check="power" ${rv.checks.power?'checked':''}> Soft-threshold behavior inspected</label><label><input type="checkbox" data-check="modules" ${rv.checks.modules?'checked':''}> ME merging and module sizes inspected</label><label><input type="checkbox" data-check="outliers" ${rv.checks.outliers?'checked':''}> Sample/outlier evidence inspected</label><label><input type="checkbox" data-check="ready" ${rv.checks.ready?'checked':''}> Ready for downstream analysis</label></div></div><div><b>Review notes</b><textarea id="reviewNotes" placeholder="Why use, hold, or reject this run?">${esc(rv.notes)}</textarea><div class="sub">Saved in this browser · ${rv.updated_at?`updated ${esc(rv.updated_at)}`:'not reviewed yet'}</div></div></section><h3 class="section-title">Visual evidence</h3><div class="artifact-tabs">${artifacts.map(a=>`<button data-artifact="${a.id}" class="${a.id===chosen?.id?'active':''}">${esc(a.label)}</button>`).join('')}</div><div class="viewer">${!chosen?'<div class="empty">No visual artifact</div>':chosen.kind==='html'?`<div style="width:100%"><div class="open"><a href="${esc(chosen.path)}" target="_blank">Open interactive result in a new tab</a></div><iframe src="${esc(chosen.path)}"></iframe></div>`:`<a href="${esc(chosen.path)}" target="_blank"><img src="${esc(chosen.path)}" alt="${esc(chosen.label)}"></a>`}</div><div class="lower"><div><h3 class="section-title">Final proper modules</h3><div class="table-wrap"><table><thead><tr><th>Module</th><th>Features</th><th>Median |kME|</th></tr></thead><tbody>${proper.map(m=>`<tr><td><span class="swatch" style="background:${esc(m.color)}"></span>${esc(m.color)}</td><td>${fmt(m.features)}</td><td>${fmt(m.median_abs_kme,3)}</td></tr>`).join('')}</tbody></table></div></div><div><h3 class="section-title">Files and provenance</h3><div class="files">${Object.entries(groups).map(([g,fs])=>`<details ${g==='summary'?'open':''}><summary>${esc(g)} (${fs.length})</summary>${fs.map(f=>`<a href="${esc(f.path)}" target="_blank">${esc(f.label)}</a>`).join('')}</details>`).join('')}<details><summary>Parameters</summary>${Object.entries(r.provenance).map(([k,v])=>`<div><b>${esc(k)}</b>: ${fmt(v)}</div>`).join('')}</details></div></div></div>`;document.querySelectorAll('[data-decision]').forEach(n=>n.onclick=()=>saveReview(r.id,{decision:n.dataset.decision}));document.querySelectorAll('[data-check]').forEach(n=>n.onchange=()=>saveReview(r.id,{checks:{...reviewOf(r.id).checks,[n.dataset.check]:n.checked}}));$('reviewNotes').onchange=e=>saveReview(r.id,{notes:e.target.value});document.querySelectorAll('[data-artifact]').forEach(n=>n.onclick=()=>{artifactChoice[r.id]=n.dataset.artifact;renderDetail()});$('addCompare').onclick=()=>{compareIds=compareIds.includes(r.id)?compareIds.filter(x=>x!==r.id):[...compareIds,r.id].slice(-4);renderAll()}}
function renderCompare(){const panel=$('comparePanel');$('compareCount').textContent=compareIds.length;if(!compareIds.length){panel.classList.remove('visible');$('compareGrid').innerHTML='';return}const selected=compareIds.map(id=>records.find(r=>r.id===id)).filter(Boolean);$('compareGrid').innerHTML=selected.map(r=>{const a=r.artifacts.find(x=>x.id==='me_dendrogram')||r.artifacts[0];return `<div class="compare-card"><h3>${esc(r.title)}</h3><div class="sub">${esc(r.dimensions.strategy)} · dynamic ${fmt(r.metrics.dynamic_modules)} → final ${fmt(r.metrics.final_modules)}</div>${a?`<a href="${esc(a.path)}" target="_blank"><img src="${esc(a.path)}" alt="${esc(a.label)}"></a>`:''}</div>`}).join('')}
function renderAll(){const view=$('viewMode').value;document.body.dataset.view=view;renderKpis();if(view==='cell'){renderCellComparison()}else{renderList();renderDetail();renderCompare()}$('sourceLine').textContent=`Schema ${bundle.schema_version} · generated ${bundle.generated_at_utc} · ${bundle.sources.map(s=>s.label).join(', ')}`}
$('showCompare').onclick=()=>{$('comparePanel').classList.toggle('visible');$('comparePanel').scrollIntoView({behavior:'smooth'})};$('clearCompare').onclick=()=>{compareIds=[];renderAll()};$('nextUnreviewed').onclick=()=>{const rows=filtered(),i=Math.max(-1,rows.findIndex(r=>r.id===selectedId));const next=[...rows.slice(i+1),...rows.slice(0,i+1)].find(r=>reviewOf(r.id).decision==='unreviewed');if(next){selectedId=next.id;renderAll()}};
function csvCell(v){const s=String(v??'');return /[",\n]/.test(s)?`"${s.replaceAll('"','""')}"`:s}function exportReview(){const head=['record_id','source_id','network_id','decision','power_checked','modules_checked','outliers_checked','ready_checked','notes','updated_at'];const rows=records.map(r=>{const v=reviewOf(r.id);return [r.id,r.source_id,r.provenance.network_id,v.decision,v.checks.power,v.checks.modules,v.checks.outliers,v.checks.ready,v.notes,v.updated_at]}).filter(x=>x[3]!=='unreviewed'||x[8]).map(x=>x.map(csvCell).join(','));const blob=new Blob([[head.join(','),...rows].join('\n')],{type:'text/csv;charset=utf-8'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=`experiment_review_${new Date().toISOString().slice(0,10)}.csv`;a.click();URL.revokeObjectURL(a.href)}$('exportReview').onclick=exportReview;$('importButton').onclick=()=>$('importReview').click();
function parseCSV(text){const rows=[];let row=[],cell='',quoted=false;for(let i=0;i<text.length;i++){const c=text[i];if(quoted&&c==='"'&&text[i+1]==='"'){cell+='"';i++}else if(c==='"'){quoted=!quoted}else if(c===','&&!quoted){row.push(cell);cell=''}else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(cell);if(row.some(x=>x!==''))rows.push(row);row=[];cell=''}else cell+=c}row.push(cell);if(row.some(x=>x!==''))rows.push(row);return rows}$('importReview').onchange=e=>{const file=e.target.files[0];if(!file)return;const reader=new FileReader();reader.onload=()=>{const rows=parseCSV(reader.result),head=rows.shift(),ix=Object.fromEntries(head.map((x,i)=>[x,i]));for(const r of rows){const id=r[ix.record_id];if(!records.some(x=>x.id===id))continue;reviews[id]={decision:r[ix.decision]||'unreviewed',checks:{power:r[ix.power_checked]==='true',modules:r[ix.modules_checked]==='true',outliers:r[ix.outliers_checked]==='true',ready:r[ix.ready_checked]==='true'},notes:r[ix.notes]||'',updated_at:r[ix.updated_at]||new Date().toISOString()}}localStorage.setItem(storageKey,JSON.stringify(reviews));renderAll()};reader.readAsText(file,'utf-8')};renderAll();
</script>
</body>
</html>'''


if __name__ == "__main__":
    arguments = parse_args()
    build(arguments.run_root.resolve())
