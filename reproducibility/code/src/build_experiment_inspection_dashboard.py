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
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import statistics
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
.final-jaccard-only { display: none; }
body[data-view="final_jaccard"] .final-jaccard-only { display: block; }
body[data-view="final_jaccard"] .cell-only,
body[data-view="final_jaccard"] .outlier-only,
body[data-view="final_jaccard"] .run-only { display: none !important; }
body[data-view="final_jaccard"] .toolbar > :not(#viewMode) { display: none !important; }
.final-jaccard-page { width: 100%; height: calc(100vh - 138px); min-height: 760px; border: 0; background: #fff; }
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
.comparison-hierarchy {
  display: grid;
  gap: 10px;
  margin: 0 0 20px;
  padding: 14px 0 16px;
  border-top: 1px solid var(--line);
  border-bottom: 1px solid var(--line);
}
.comparison-dimension {
  display: grid;
  grid-template-columns: 118px minmax(0, 1fr);
  gap: 16px;
  align-items: center;
}
.comparison-dimension-label {
  color: var(--text-muted);
  font-size: .74rem;
  font-weight: 600;
  letter-spacing: .04em;
  text-transform: uppercase;
}
.comparison-dimension-options {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
}
.comparison-dimension-options button {
  min-height: 34px;
  padding: 6px 13px;
  border: 1px solid var(--line-strong);
  background: var(--surface);
  color: var(--text-muted);
  font: inherit;
  font-size: .8rem;
  cursor: pointer;
}
.comparison-dimension-options button:hover { border-color: var(--text-muted); color: var(--text); }
.comparison-dimension-options button.active {
  border-color: var(--text);
  background: var(--text);
  color: var(--surface);
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
.strategy-grid { grid-template-columns: repeat(3, minmax(0, 1fr)); }
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
  .strategy-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
}
@media (max-width: 720px) {
  .comparison-dimension { grid-template-columns: 1fr; gap: 6px; }
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
  .strategy-grid { grid-template-columns: 1fr; }
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


OUTLIER_REVIEW_CSS = r"""
body[data-view="outliers"] .cell-only,
body[data-view="outliers"] .run-only { display: none !important; }
body:not([data-view="outliers"]) .outlier-only { display: none !important; }
body[data-view="outliers"] .toolbar > :not(#viewMode) { display: none !important; }
.outlier-page { max-width: var(--page-max); margin: 0 auto; padding: 31px 32px 64px; border-left: 1px solid var(--line-subtle); border-right: 1px solid var(--line-subtle); }
.outlier-page .comparison-level { display: grid; grid-template-columns: 190px minmax(0, 1fr); gap: 16px; align-items: start; }
.comparison-level-label { display: flex; align-items: flex-start; gap: 9px; }
.comparison-level-label > span { display: grid; place-items: center; width: 22px; height: 22px; border: 1px solid var(--line-strong); border-radius: 50%; color: var(--text-muted); font-size: .7rem; }
.comparison-level-label b { display: block; font-size: .78rem; font-weight: 600; }
.comparison-level-label small { display: block; margin-top: 2px; color: var(--text-muted); font-size: .68rem; }
.comparison-switches { display: flex; flex-wrap: wrap; gap: 6px; }
.comparison-switches button { min-width: 116px; padding: 7px 11px; border: 1px solid var(--line-strong); background: var(--surface); color: var(--text-muted); text-align: left; }
.comparison-switches button span { display: block; color: inherit; font-size: .78rem; }
.comparison-switches button small { display: block; margin-top: 2px; color: inherit; opacity: .75; font-size: .66rem; }
.comparison-switches button.active { border-color: var(--text); background: var(--text); color: var(--surface); }
.outlier-context, .outlier-method-note, .outlier-matrix-note { margin: 10px 0 16px; padding: 10px 12px; border-left: 3px solid var(--line-strong); background: var(--surface-subtle); color: var(--text-muted); font-size: .76rem; line-height: 1.55; }
.outlier-context button { margin-left: 12px; padding: 4px 8px; }
.outlier-controls { display: grid; grid-template-columns: minmax(260px, 1fr) minmax(160px, .45fr) minmax(180px, .5fr); gap: 10px; margin: 18px 0 10px; }
.outlier-controls label span { display: block; margin-bottom: 4px; color: var(--text-muted); font-size: .7rem; }
.outlier-controls input, .outlier-controls select { width: 100%; padding: 8px 9px; }
.outlier-downloads { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: center; margin-bottom: 12px; font-size: .74rem; }
.outlier-downloads a { color: var(--link); }
.outlier-downloads button { margin-left: auto; padding: 6px 9px; }
.outlier-candidate-table { max-height: 560px; }
.outlier-candidate-table th, .outlier-candidate-table td { padding: 7px 9px; text-align: center; }
.outlier-candidate-table th:nth-child(2), .outlier-candidate-table td:nth-child(2) { position: sticky; left: 0; z-index: 2; min-width: 235px; background: var(--surface); text-align: left; }
.outlier-candidate-table thead th:nth-child(2) { z-index: 3; }
.outlier-candidate-table tr.selected td { background: var(--selected); }
.outlier-candidate-table tr { cursor: pointer; }
.outlier-cell-head, .outlier-cell { min-width: 92px; }
.outlier-cell-head { white-space: normal; }
.outlier-cell.flagged { background: #fff1f1; color: #991b1b; }
.outlier-cell.normal { background: #f1f8f4; color: #17633a; }
.outlier-cell.na { background: var(--surface-subtle); color: var(--text-muted); }
.outlier-cell-status { display: flex; justify-content: center; gap: 5px; }
.outlier-cell small { display: block; margin-top: 2px; color: inherit; opacity: .8; font-size: .64rem; }
.table-secondary { margin-top: 2px; color: var(--text-muted); font-size: .67rem; }
.outlier-decision { display: inline-block; padding: 2px 6px; border: 1px solid var(--line-strong); color: var(--text-muted); font-size: .66rem; text-transform: uppercase; }
.outlier-decision.keep { color: var(--ok); border-color: var(--ok); }
.outlier-decision.hold { color: var(--warn); border-color: var(--warn); }
.outlier-decision.exclude { color: var(--danger); border-color: var(--danger); }
.outlier-detail { margin-top: 28px; padding-top: 18px; border-top: 1px solid var(--line-strong); }
.outlier-detail-head { display: flex; align-items: flex-start; justify-content: space-between; gap: 20px; }
.outlier-detail-head h2 { margin: 0 0 3px; font-size: 1.08rem; font-weight: 500; }
.outlier-evidence-grid { display: grid; grid-template-columns: repeat(4, 1fr); margin: 16px 0; border-top: 1px solid var(--line); border-bottom: 1px solid var(--line); }
.outlier-evidence-grid > div { padding: 10px 12px; border-right: 1px solid var(--line-subtle); }
.outlier-evidence-grid > div:last-child { border-right: 0; }
.outlier-evidence-grid span { display: block; color: var(--text-muted); font-size: .68rem; }
.outlier-evidence-grid b { display: block; margin-top: 2px; font-size: 1rem; font-weight: 500; }
.outlier-decision-panel { display: grid; grid-template-columns: minmax(300px, .8fr) minmax(320px, 1.2fr); gap: 20px; margin: 14px 0; padding: 14px; background: var(--surface-subtle); }
.outlier-decision-panel textarea { width: 100%; min-height: 72px; padding: 8px; resize: vertical; }
.outlier-cell-evidence { max-height: 420px; }
.outlier-cell-inline { display: inline-block; min-width: 26px; padding: 2px 6px; text-align: center; }
.outlier-cell-inline.flagged { background: #fff1f1; color: #991b1b; }
.outlier-cell-inline.normal { background: #f1f8f4; color: #17633a; }
.outlier-cell-inline.na { background: var(--surface-subtle); color: var(--text-muted); }
.outlier-link { padding: 4px 7px; }
@media (max-width: 720px) {
  .outlier-page .comparison-level { grid-template-columns: 1fr; }
  .outlier-controls, .outlier-decision-panel { grid-template-columns: 1fr; }
  .outlier-evidence-grid { grid-template-columns: repeat(2, 1fr); }
}
@media (max-width: 820px) { .outlier-page { padding: 25px 18px 48px; border-left: 0; border-right: 0; } }
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
        after_module_root = (
            network_root / "outlier_removal/05_find_modules" / MODULE_CONDITION
        )
        after_analysis_path = after_module_root / "analysis_summary.json"
        after_analysis = read_json(after_analysis_path) if after_analysis_path.exists() else {}
        after_diagnosis_root = (
            network_root / "outlier_removal/06_kme_module_diagnosis" / MODULE_CONDITION
        )
        after_diagnosis_summary_path = after_diagnosis_root / "analysis_summary.json"
        after_diagnosis_summary = (
            read_json(after_diagnosis_summary_path)
            if after_diagnosis_summary_path.exists()
            else {}
        )
        after_pairwise_summary_path = (
            after_diagnosis_root / "pairwise_hd_analysis_summary.json"
        )
        after_pairwise_summary = (
            read_json(after_pairwise_summary_path)
            if after_pairwise_summary_path.exists()
            else {}
        )
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
        after_cross_network_root = (
            run_root
            / "experiments/05_cross_omics_module_similarity/after_outlier_removal/transcriptomics"
            / row["threshold"]
            / row["strategy"]
            / row["mode"]
            / row["specimen"]
        )
        if row["mode"] == "cell_separate":
            after_cross_network_root = after_cross_network_root / row["cell_type"]
        after_cross_summary_path = after_cross_network_root / "analysis_summary.json"
        after_cross_summary = (
            read_json(after_cross_summary_path) if after_cross_summary_path.exists() else {}
        )
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
        mirrored_image_artifact(
            artifacts,
            after_diagnosis_root / "figures/module_diagnosis_pairwise_hd_heatmap.png",
            dashboard_root,
            "after_outlier_diagnosis_pairwise_hd",
            "Module-diagnosis (vs HD)",
            "module diagnosis after outlier removal",
        )
        for spec in (
            ("me_dendrogram", "ME dendrogram", "05_find_modules/" + MODULE_CONDITION + "/figures/module_eigengene_dendrogram.png", "findModules", "image"),
            ("gene_dendrogram", "Gene dendrogram", "05_find_modules/" + MODULE_CONDITION + "/figures/gene_dendrogram_modules.png", "findModules", "image"),
            ("soft_threshold", "Soft threshold", "05_find_modules/" + MODULE_CONDITION + "/figures/soft_threshold.png", "findModules", "image"),
            ("module_sizes", "Module sizes", "05_find_modules/" + MODULE_CONDITION + "/figures/module_sizes.png", "findModules", "image"),
            ("diagnosis_heatmap", "One-hot diagnosis (audit)", "06_kme_module_diagnosis/" + MODULE_CONDITION + "/figures/module_diagnosis_heatmap.png", "module diagnosis", "image"),
            ("feature_selection", "Feature selection", "01_feature_selection/feature_selection.png", "feature selection", "image"),
            ("before_heatmap", "Heatmap", "02_before_preprocess/02_data_distribution/heatmap.png", "before preprocess", "image"),
            ("before_pca", "PCA metadata", "02_before_preprocess/02_data_distribution/pca_metadata.png", "before preprocess", "image"),
            ("before_pca_3d", "PCA 3D interactive", "02_before_preprocess/02_data_distribution/pca_3d_interactive.html", "before preprocess", "html"),
            ("after_outlier_heatmap", "Heatmap", "outlier_removal/02_before_preprocess/02_data_distribution/heatmap.png", "before preprocess", "image"),
            ("after_outlier_pca", "PCA metadata", "outlier_removal/02_before_preprocess/02_data_distribution/pca_metadata.png", "before preprocess", "image"),
            ("after_outlier_pca_3d", "PCA 3D interactive", "outlier_removal/02_before_preprocess/02_data_distribution/pca_3d_interactive.html", "before preprocess", "html"),
            ("cut_diagnostic", "Cut diagnostic", "03_preprocess/01_cut_diagnostic/dendrogram_cut_diagnostic.png", "preprocess", "image"),
            ("selected_cut", "Selected cut", "03_preprocess/02_selected_cut/dendrogram_selected_cut.png", "preprocess", "image"),
            ("after_outlier_cut_diagnostic", "Cut diagnostic", "outlier_removal/03_preprocess/01_cut_diagnostic/dendrogram_cut_diagnostic.png", "preprocess", "image"),
            ("after_outlier_me_dendrogram", "ME dendrogram", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/figures/module_eigengene_dendrogram.png", "findModules after outlier removal", "image"),
            ("after_outlier_gene_dendrogram", "Gene dendrogram", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/figures/gene_dendrogram_modules.png", "findModules after outlier removal", "image"),
            ("after_outlier_soft_threshold", "Soft threshold", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/figures/soft_threshold.png", "findModules after outlier removal", "image"),
            ("after_outlier_module_sizes", "Module sizes", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/figures/module_sizes.png", "findModules after outlier removal", "image"),
            ("after_outlier_diagnosis_heatmap", "One-hot diagnosis (audit)", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/figures/module_diagnosis_heatmap.png", "module diagnosis after outlier removal", "image"),
            ("after_heatmap", "After heatmap", "04_after_preprocess/02_data_distribution/heatmap.png", "after preprocess", "image"),
            ("after_pca", "After PCA", "04_after_preprocess/02_data_distribution/pca_metadata.png", "after preprocess", "image"),
            ("after_pca_3d", "After PCA 3D", "04_after_preprocess/02_data_distribution/pca_3d_interactive.html", "after preprocess", "html"),
            ("after_outlier_after_heatmap", "After heatmap", "outlier_removal/04_after_preprocess/02_data_distribution/heatmap.png", "after preprocess after outlier removal", "image"),
            ("after_outlier_after_pca", "After PCA", "outlier_removal/04_after_preprocess/02_data_distribution/pca_metadata.png", "after preprocess after outlier removal", "image"),
            ("after_outlier_after_pca_3d", "After PCA 3D", "outlier_removal/04_after_preprocess/02_data_distribution/pca_3d_interactive.html", "after preprocess after outlier removal", "html"),
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
        existing_artifact(
            artifacts,
            after_cross_network_root,
            dashboard_root,
            "after_outlier_cross_omics_similarity",
            "RNA-protein module similarity",
            "kme_weighted_jaccard_heatmap.png",
            "cross-omics after outlier removal",
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
            ("outlier_removal_excluded_samples", "Confirmed excluded samples", "outlier_removal/01_sample_selection/excluded_samples.csv", "outlier removal"),
            ("outlier_removal_selection_summary", "Outlier-removal selection summary", "outlier_removal/01_sample_selection/selection_summary.json", "outlier removal"),
            ("after_outlier_pca_scores", "After-outlier-removal PCA scores", "outlier_removal/02_before_preprocess/02_data_distribution/pca_scores.csv", "outlier removal QC"),
            ("after_outlier_after_pca_scores", "After-outlier post-preprocess PCA scores", "outlier_removal/04_after_preprocess/02_data_distribution/pca_scores.csv", "after preprocess after outlier removal"),
            ("after_outlier_cut_recommendation", "Largest-gap cut recommendation", "outlier_removal/03_preprocess/01_cut_diagnostic/cut_recommendation.json", "outlier removal preprocess"),
            ("after_outlier_cut_candidates", "All cut candidates", "outlier_removal/03_preprocess/01_cut_diagnostic/cut_candidates.csv", "outlier removal preprocess"),
            ("after_outlier_preprocess_summary", "cut=inf preprocess summary", "outlier_removal/03_preprocess/01_cut_diagnostic/preprocess_summary.json", "outlier removal preprocess"),
            ("after_outlier_run_manifest", "Outlier-removal run manifest", "outlier_removal/run_manifest.json", "provenance"),
            ("after_outlier_analysis_summary", "findModules analysis summary", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/analysis_summary.json", "summary after outlier removal"),
            ("after_outlier_module_summary", "Module summary", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/module_summary.csv", "module after outlier removal"),
            ("after_outlier_module_assignments", "Module assignments and kME", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/module_assignments.csv.gz", "module after outlier removal"),
            ("after_outlier_dynamic_mes", "Dynamic module eigengenes", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/dynamic_module_eigengenes.csv", "module after outlier removal"),
            ("after_outlier_final_mes", "Final module eigengenes", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/module_eigengenes.csv", "module after outlier removal"),
            ("after_outlier_soft_threshold_table", "Soft-threshold table", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/soft_threshold_table.csv", "module after outlier removal"),
            ("after_outlier_findmodules_manifest", "findModules run manifest", "outlier_removal/05_find_modules/" + MODULE_CONDITION + "/run_manifest.json", "provenance"),
            ("after_outlier_official_kme", "Official all-module kME", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/feature_module_kme.csv.gz", "kME after outlier removal"),
            ("after_outlier_kme_summary", "Official kME summary", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_kme_summary.csv", "kME after outlier removal"),
            ("after_outlier_diagnosis_statistics", "One-hot diagnosis statistics (audit)", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_statistics.csv", "module diagnosis after outlier removal"),
            ("after_outlier_diagnosis_pairwise_statistics", "Pairwise disease-vs-HD statistics", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/module_diagnosis_pairwise_hd_statistics.csv", "module diagnosis after outlier removal"),
            ("after_outlier_diagnosis_pairwise_summary", "Pairwise disease-vs-HD analysis summary", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/pairwise_hd_analysis_summary.json", "summary after outlier removal"),
            ("after_outlier_diagnosis_manifest", "Official kME/diagnosis run manifest", "outlier_removal/06_kme_module_diagnosis/" + MODULE_CONDITION + "/run_manifest.json", "provenance"),
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
            ("cross_all_gene_contributions", "All module-pair overlapping genes", "all_pair_gene_contributions.csv.gz", "cross-omics"),
            ("cross_gene_contributions", "Reciprocal-best gene contributions", "reciprocal_best_gene_contributions.csv.gz", "cross-omics"),
            ("cross_mapping_audit", "Cross-omics gene mapping audit", "gene_mapping_audit.json", "cross-omics"),
            ("cross_summary", "Cross-omics analysis summary", "analysis_summary.json", "cross-omics"),
            ("cross_manifest", "Cross-omics run manifest", "run_manifest.json", "provenance"),
        ):
            existing_file(files, cross_network_root, dashboard_root, *spec)
        for spec in (
            ("after_outlier_cross_pair_similarity", "All RNA-protein module pairs", "module_pair_similarity.csv", "cross-omics after outlier removal"),
            ("after_outlier_cross_best_matches", "Transcript module best matches", "transcript_module_best_matches.csv", "cross-omics after outlier removal"),
            ("after_outlier_cross_all_gene_contributions", "All module-pair overlapping genes", "all_pair_gene_contributions.csv.gz", "cross-omics after outlier removal"),
            ("after_outlier_cross_gene_contributions", "Reciprocal-best gene contributions", "reciprocal_best_gene_contributions.csv.gz", "cross-omics after outlier removal"),
            ("after_outlier_cross_mapping_audit", "Cross-omics gene mapping audit", "gene_mapping_audit.json", "cross-omics after outlier removal"),
            ("after_outlier_cross_summary", "Cross-omics analysis summary", "analysis_summary.json", "cross-omics after outlier removal"),
            ("after_outlier_cross_manifest", "Cross-omics run manifest", "run_manifest.json", "provenance"),
        ):
            existing_file(files, after_cross_network_root, dashboard_root, *spec)

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
                    "after_outlier_samples": after_analysis.get("samples"),
                    "after_outlier_genes": after_analysis.get("genes"),
                    "after_outlier_selected_power": after_analysis.get("selected_power"),
                    "after_outlier_dynamic_modules": after_analysis.get("dynamic_module_count"),
                    "after_outlier_final_modules": after_analysis.get("final_module_count"),
                    "after_outlier_unassigned_fraction": after_analysis.get("raw_unassigned_fraction"),
                    "after_outlier_largest_module_fraction": after_analysis.get("largest_module_fraction"),
                    "after_outlier_median_abs_kme": after_analysis.get("median_assigned_abs_kME"),
                    "after_outlier_runtime_seconds": after_analysis.get("runtime_seconds"),
                    "after_outlier_finite_cut_applied": after_analysis.get("finite_sample_cut_applied"),
                    "after_outlier_cross_omics_common_genes": after_cross_summary.get("common_gene_universe"),
                    "after_outlier_cross_omics_max_similarity": after_cross_summary.get("maximum_weighted_jaccard"),
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
                    "after_outlier_findmodules_status": after_analysis.get("status"),
                    "after_outlier_kme_method": after_diagnosis_summary.get("kme_method"),
                    "after_outlier_module_diagnosis_method": after_pairwise_summary.get("method") or after_diagnosis_summary.get("heatmap_method"),
                    "after_outlier_cross_omics_method": "absolute own-module kME weighted Jaccard" if after_cross_summary else None,
                    "after_outlier_sample_cut": "none; cut=inf preprocess matrix used directly",
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


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _open_csv(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def build_module_member_shards(
    records: list[dict[str, Any]], dashboard_root: Path
) -> dict[str, Any]:
    """Build lazy, file:// compatible module-member data for the explorer.

    Each module is a separate JavaScript shard.  The browser therefore loads
    only the selected module rather than every feature across every network.
    Values come directly from the official CalculateSignedKME output; missing
    kME values remain JSON null instead of being converted to zero.
    """

    member_root = dashboard_root / "module-members"
    shard_root = member_root / "shards"
    shard_root.mkdir(parents=True, exist_ok=True)
    index: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "records": {},
    }
    expected_shards: set[Path] = set()
    total_rows = 0

    def add_member_record(
        record: dict[str, Any],
        kme_file: dict[str, Any],
        member_record_id: str,
        member_title: str,
    ) -> None:
        nonlocal total_rows
        source_path = (dashboard_root / kme_file["path"]).resolve()
        if not source_path.exists():
            return

        grouped: dict[str, list[list[Any]]] = {}
        with _open_csv(source_path) as handle:
            for row in csv.DictReader(handle):
                module_color = (
                    row.get("module_color")
                    or row.get("dynamic_color")
                    or "unassigned"
                ).strip()
                raw_unassigned = _truthy(row.get("is_raw_unassigned"))
                signed_kme = parse_number(row.get("official_assigned_signed_kME"))
                if signed_kme is None and not raw_unassigned:
                    signed_kme = parse_number(row.get("official_module_kME"))
                abs_kme = (
                    None
                    if raw_unassigned
                    else parse_number(row.get("official_assigned_abs_kME"))
                )
                if abs_kme is None and signed_kme is not None:
                    abs_kme = abs(float(signed_kme))
                rank = None if raw_unassigned else parse_number(row.get("hub_rank_within_module"))
                grouped.setdefault(module_color, []).append(
                    [
                        row.get("feature_id") or row.get("gene_id") or row.get("PROBEID") or "",
                        row.get("gene_id") or "",
                        row.get("SYMBOL") or "",
                        row.get("PROBEID") or "",
                        row.get("UNIPROT") or "",
                        signed_kme,
                        abs_kme,
                        rank,
                        raw_unassigned,
                    ]
                )

        module_index: dict[str, Any] = {}
        for module_color, rows in sorted(grouped.items()):
            rows.sort(
                key=lambda row: (
                    row[7] is None,
                    row[7] if row[7] is not None else float("inf"),
                    -(row[6] if row[6] is not None else -1),
                    str(row[0]),
                )
            )
            key = f"{member_record_id}::{module_color}"
            filename = hashlib.sha256(key.encode("utf-8")).hexdigest()[:24] + ".js"
            shard_path = shard_root / filename
            payload = {
                "columns": [
                    "feature_id",
                    "gene_id",
                    "symbol",
                    "probe_id",
                    "uniprot",
                    "signed_kme",
                    "abs_kme",
                    "hub_rank",
                    "is_raw_unassigned",
                ],
                "rows": rows,
            }
            shard_path.write_text(
                "window.MODULE_MEMBER_SHARDS=window.MODULE_MEMBER_SHARDS||{};"
                f"window.MODULE_MEMBER_SHARDS[{json.dumps(key, ensure_ascii=False)}]="
                + json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
                + ";\n",
                encoding="utf-8",
            )
            expected_shards.add(shard_path.resolve())
            abs_values = [float(row[6]) for row in rows if row[6] is not None]
            raw_unassigned = sum(bool(row[8]) for row in rows)
            module_index[module_color] = {
                "path": f"module-members/shards/{filename}",
                "count": len(rows),
                "assigned_kme_count": len(abs_values),
                "raw_unassigned_count": raw_unassigned,
                "fully_raw_unassigned": raw_unassigned == len(rows),
                "median_abs_kme": statistics.median(abs_values) if abs_values else None,
                "max_abs_kme": max(abs_values) if abs_values else None,
            }
            total_rows += len(rows)

        index["records"][member_record_id] = {
            "title": member_title,
            "source_id": record.get("source_id"),
            "network_id": record.get("provenance", {}).get("network_id"),
            "source_path": kme_file["path"],
            "total_features": sum(item["count"] for item in module_index.values()),
            "modules": module_index,
        }

    for record in records:
        files_by_id = {item.get("id"): item for item in record.get("files", [])}
        official_kme = files_by_id.get("official_kme")
        if official_kme:
            add_member_record(
                record,
                official_kme,
                record["id"],
                record.get("title") or record["id"],
            )
        after_outlier_kme = files_by_id.get("after_outlier_official_kme")
        if after_outlier_kme:
            add_member_record(
                record,
                after_outlier_kme,
                f"{record['id']}::after_outlier",
                f"{record.get('title', record['id'])} - After outlier removal",
            )

    for stale in shard_root.glob("*.js"):
        if stale.resolve() not in expected_shards:
            stale.unlink()

    index_path = member_root / "module-member-index.js"
    index_path.write_text(
        "window.MODULE_MEMBER_INDEX="
        + json.dumps(index, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )
    return {
        "records": len(index["records"]),
        "modules": len(expected_shards),
        "features": total_rows,
        "index_path": index_path,
    }


def build_cross_omics_overlap_shards(
    records: list[dict[str, Any]], dashboard_root: Path
) -> dict[str, Any]:
    """Build lazy offline data for inspecting genes shared by each module pair."""

    overlap_root = dashboard_root / "cross-omics-overlap"
    shard_root = overlap_root / "shards"
    shard_root.mkdir(parents=True, exist_ok=True)
    index: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "records": {},
    }
    expected_shards: set[Path] = set()
    total_pairs = 0
    total_rows = 0

    overlap_records: list[dict[str, Any]] = []
    for record in records:
        if str(record.get("dimensions", {}).get("cell_type") or "").strip().lower() == "other":
            continue
        overlap_records.append(record)
        files_by_id = {item.get("id"): item for item in record.get("files", [])}
        after_pair = files_by_id.get("after_outlier_cross_pair_similarity")
        after_contributions = files_by_id.get(
            "after_outlier_cross_all_gene_contributions"
        )
        if after_pair and after_contributions:
            after_record = {
                **record,
                "id": f"{record['id']}::after_outlier",
                "title": f"{record.get('title', record['id'])} - After outlier removal",
                "files": [
                    {**after_pair, "id": "cross_pair_similarity"},
                    {
                        **after_contributions,
                        "id": "cross_all_gene_contributions",
                    },
                ],
                "provenance": {
                    **record.get("provenance", {}),
                    "outlier_removal_state": "after",
                },
            }
            overlap_records.append(after_record)

    for record in overlap_records:
        files = {item.get("id"): item for item in record.get("files", [])}
        pair_file = files.get("cross_pair_similarity")
        contribution_file = files.get("cross_all_gene_contributions")
        if not pair_file or not contribution_file:
            continue
        pair_path = (dashboard_root / pair_file["path"]).resolve()
        contribution_path = (dashboard_root / contribution_file["path"]).resolve()
        if not pair_path.is_file() or not contribution_path.is_file():
            continue

        pairs: list[dict[str, Any]] = []
        with _open_csv(pair_path) as handle:
            for row in csv.DictReader(handle):
                pairs.append(
                    {
                        "transcript_module": row.get("transcript_module") or "",
                        "proteomics_module": row.get("proteomics_module") or "",
                        "transcript_module_genes": parse_number(
                            row.get("transcript_module_genes_in_common")
                        ),
                        "proteomics_module_genes": parse_number(
                            row.get("proteomics_module_genes_in_common")
                        ),
                        "shared_genes": parse_number(row.get("shared_genes")),
                        "jaccard": parse_number(row.get("jaccard")),
                        "weighted_intersection": parse_number(
                            row.get("weighted_intersection")
                        ),
                        "weighted_union": parse_number(row.get("weighted_union")),
                        "weighted_jaccard": parse_number(
                            row.get("kme_weighted_jaccard")
                        ),
                        "transcript_best": _truthy(row.get("transcript_best")),
                        "proteomics_best": _truthy(row.get("proteomics_best")),
                        "reciprocal_best": _truthy(row.get("reciprocal_best")),
                    }
                )
        pairs.sort(
            key=lambda item: (
                -(float(item["weighted_jaccard"]) if item["weighted_jaccard"] is not None else -1),
                item["transcript_module"],
                item["proteomics_module"],
            )
        )

        rows: list[list[Any]] = []
        with _open_csv(contribution_path) as handle:
            for row in csv.DictReader(handle):
                rows.append(
                    [
                        row.get("transcript_module") or "",
                        row.get("proteomics_module") or "",
                        row.get("gene") or "",
                        parse_number(row.get("transcript_abs_kME")),
                        parse_number(row.get("proteomics_abs_kME")),
                        parse_number(row.get("weighted_intersection")),
                        parse_number(row.get("weighted_union")),
                    ]
                )
        rows.sort(
            key=lambda row: (
                str(row[0]),
                str(row[1]),
                -(float(row[5]) if row[5] is not None else -1),
                str(row[2]),
            )
        )

        filename = hashlib.sha256(record["id"].encode("utf-8")).hexdigest()[:24] + ".js"
        shard_path = shard_root / filename
        payload = {
            "columns": [
                "transcript_module",
                "proteomics_module",
                "gene",
                "transcript_abs_kme",
                "proteomics_abs_kme",
                "weighted_intersection",
                "weighted_union",
            ],
            "rows": rows,
        }
        shard_path.write_text(
            "window.CROSS_OMICS_OVERLAP_SHARDS=window.CROSS_OMICS_OVERLAP_SHARDS||{};"
            f"window.CROSS_OMICS_OVERLAP_SHARDS[{json.dumps(record['id'], ensure_ascii=False)}]="
            + json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
            + ";\n",
            encoding="utf-8",
        )
        expected_shards.add(shard_path.resolve())
        transcript_modules = list(
            dict.fromkeys(item["transcript_module"] for item in pairs)
        )
        proteomics_modules = list(
            dict.fromkeys(item["proteomics_module"] for item in pairs)
        )
        index["records"][record["id"]] = {
            "title": record.get("title"),
            "network_id": record.get("provenance", {}).get("network_id"),
            "specimen": record.get("dimensions", {}).get("specimen"),
            "cell_type": record.get("dimensions", {}).get("cell_type"),
            "threshold": record.get("dimensions", {}).get("threshold"),
            "strategy": record.get("dimensions", {}).get("strategy"),
            "mode": record.get("dimensions", {}).get("mode"),
            "proteomics_condition": record.get("provenance", {}).get(
                "cross_omics_proteomics_condition"
            ),
            "path": f"cross-omics-overlap/shards/{filename}",
            "transcript_modules": transcript_modules,
            "proteomics_modules": proteomics_modules,
            "pair_columns": [
                "transcript_module",
                "proteomics_module",
                "transcript_module_genes",
                "proteomics_module_genes",
                "shared_genes",
                "jaccard",
                "weighted_intersection",
                "weighted_union",
                "weighted_jaccard",
                "transcript_best",
                "proteomics_best",
                "reciprocal_best",
            ],
            "pairs": [
                [
                    item["transcript_module"],
                    item["proteomics_module"],
                    item["transcript_module_genes"],
                    item["proteomics_module_genes"],
                    item["shared_genes"],
                    item["jaccard"],
                    item["weighted_intersection"],
                    item["weighted_union"],
                    item["weighted_jaccard"],
                    item["transcript_best"],
                    item["proteomics_best"],
                    item["reciprocal_best"],
                ]
                for item in pairs
            ],
            "overlap_rows": len(rows),
        }
        total_pairs += len(pairs)
        total_rows += len(rows)

    for stale in shard_root.glob("*.js"):
        if stale.resolve() not in expected_shards:
            stale.unlink()

    index_path = overlap_root / "overlap-index.js"
    index_path.write_text(
        "window.CROSS_OMICS_OVERLAP_INDEX="
        + json.dumps(index, ensure_ascii=False, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )
    return {
        "records": len(index["records"]),
        "pairs": total_pairs,
        "overlap_rows": total_rows,
        "shards": len(expected_shards),
        "index_path": index_path,
    }


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
    module_members = build_module_member_shards(records, dashboard_root)
    cross_omics_overlap = build_cross_omics_overlap_shards(records, dashboard_root)
    bundle = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "run_root": run_root.as_posix(),
        "sources": sources,
        "records": records,
        "module_member_explorer": {
            "records": module_members["records"],
            "modules": module_members["modules"],
            "features": module_members["features"],
            "index": "module-members/module-member-index.js",
        },
        "cross_omics_overlap_explorer": {
            "records": cross_omics_overlap["records"],
            "pairs": cross_omics_overlap["pairs"],
            "overlap_rows": cross_omics_overlap["overlap_rows"],
            "index": "cross-omics-overlap/overlap-index.js",
        },
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
        DESIGN_SYSTEM_CSS
        + "\n"
        + OUTLIER_REVIEW_CSS
        + "\n"
        + MODULE_EXPLORER_CSS
        + "\n"
        + CROSS_OMICS_EXPLORER_CSS
        + "\n</style>",
        1,
    ).replace(
        '<script src="dashboard-data.js"></script>',
        MODULE_EXPLORER_HTML
        + CROSS_OMICS_EXPLORER_HTML
        + '\n<script src="module-members/module-member-index.js"></script>'
        + '\n<script src="module-members/cross-omics-index.js"></script>'
        + '\n<script src="module-members/cross-omics-view.js"></script>'
        + '\n<script src="cross-omics-overlap/overlap-index.js"></script>'
        + '\n<script src="outlier-candidates/outlier-candidates.js"></script>'
        + '\n<script src="outlier-candidates/outlier-review.js"></script>'
        + '\n<script src="dashboard-data.js"></script>',
        1,
    ).replace(
        "</body>",
        "<script>\n"
        + MODULE_EXPLORER_JS
        + "\n"
        + CROSS_OMICS_EXPLORER_JS
        + "\n</script>\n</body>",
        1,
    )
    html_path.write_text(dashboard_html, encoding="utf-8")
    shutil.copyfile(
        Path(__file__).with_name("final_jaccard_dashboard.html"),
        dashboard_root / "final-jaccard.html",
    )
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
    explorer_readme_path = dashboard_root / "MODULE_EXPLORER.md"
    explorer_readme_path.write_text(
        f"""# Module Explorer

`index.html`의 `Modules & kME` 버튼 또는 Individual run의 module table 행을 누르면
우측 Module Explorer가 열린다. 별도 web server 없이 `file://`에서 동작한다.

- indexed records: {module_members['records']}
- module shards: {module_members['modules']}
- indexed features: {module_members['features']}
- source: 각 record의 공식 `feature_module_kme.csv.gz`
- default order: `hub_rank_within_module` 오름차순
- shown values: signed kME, absolute kME, hub rank, feature/gene/protein identifiers
- raw-unassigned: 목록에는 유지하되 kME와 rank는 `NA`

검색, 정렬, 페이지 이동, 현재 검색 결과 CSV export를 지원한다. kME는 module membership
지표이며 통계적 유의성이나 질병 연관성의 p-value가 아니다.

재생성:

```powershell
python results_20261119/code/src/build_experiment_inspection_dashboard.py --run-root results_20261119
```
""",
        encoding="utf-8",
    )
    cross_overlap_readme_path = dashboard_root / "CROSS_OMICS_GENE_OVERLAP.md"
    cross_overlap_readme_path.write_text(
        f"""# RNA-protein Module Comparison Explorer

Jaccard heatmap의 transcriptomics module과 proteomics module pair를 선택해 실제 공통
유전자만 확인하는 offline explorer이다.

- indexed transcriptomics networks: {cross_omics_overlap['records']}
- indexed module pairs: {cross_omics_overlap['pairs']}
- overlapping gene rows: {cross_omics_overlap['overlap_rows']}
- source: `all_pair_gene_contributions.csv.gz`
- shown values: gene, transcriptomics |kME|, proteomics |kME|,
  `min(|kME|)`, `max(|kME|)`
- controls: module-pair selection, gene search, sorting, pagination, CSV export

표에는 두 module에 모두 존재하는 유전자만 나타난다. Weighted Jaccard의 denominator는
한쪽 module에만 존재하는 유전자도 포함하므로 표의 `max(|kME|)` 합이 전체 weighted
union과 같다는 의미는 아니다. 이 similarity는 descriptive overlap이며 유의성 검정이 아니다.
""",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": generated_at,
        "records": len(records),
        "sources": sources,
        "outputs": [
            html_path.name,
            data_path.name,
            readme_path.name,
            explorer_readme_path.name,
            cross_overlap_readme_path.name,
            module_members["index_path"].relative_to(dashboard_root).as_posix(),
            cross_omics_overlap["index_path"].relative_to(dashboard_root).as_posix(),
        ],
        "module_member_explorer": {
            "records": module_members["records"],
            "modules": module_members["modules"],
            "features": module_members["features"],
        },
        "cross_omics_overlap_explorer": {
            "records": cross_omics_overlap["records"],
            "pairs": cross_omics_overlap["pairs"],
            "overlap_rows": cross_omics_overlap["overlap_rows"],
            "shards": cross_omics_overlap["shards"],
        },
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


MODULE_EXPLORER_CSS = r"""
.module-explorer-overlay[hidden] { display: none; }
.module-explorer-overlay {
  position: fixed; inset: 0; z-index: 1000;
  background: rgba(20, 20, 20, .38);
  display: flex; justify-content: flex-end;
}
.module-explorer {
  width: min(1040px, 96vw); height: 100vh; overflow: auto;
  background: var(--surface); border-left: 1px solid var(--line-strong);
  box-shadow: -12px 0 36px rgba(0,0,0,.12); padding: 24px 28px 40px;
}
.module-explorer-head { display: flex; gap: 20px; align-items: flex-start; border-bottom: 1px solid var(--line); padding-bottom: 15px; }
.module-explorer-head h2 { margin: 0 0 3px; font-size: 23px; }
.module-explorer-head button { margin-left: auto; padding: 7px 12px; }
.module-explorer-controls { display: grid; grid-template-columns: minmax(190px,.8fr) minmax(220px,1fr) minmax(160px,.7fr) auto; gap: 8px; margin: 16px 0 10px; }
.module-explorer-controls select, .module-explorer-controls input, .module-explorer-controls button { min-width: 0; padding: 8px 9px; }
.module-explorer-summary { display: flex; flex-wrap: wrap; gap: 7px 16px; padding: 10px 12px; background: var(--panel); border: 1px solid var(--line); font-size: 13px; }
.module-explorer-note { color: var(--text-muted); font-size: 12px; margin: 9px 0 13px; }
.module-member-table { max-height: calc(100vh - 330px); }
.module-member-table tr.raw-unassigned { background: #faf7f1; }
.module-member-table td.num { text-align: right; font-variant-numeric: tabular-nums; }
.module-member-table button { border: 0; border-bottom: 1px solid var(--line-strong); background: transparent; padding: 0; }
.raw-badge { display: inline-block; padding: 1px 5px; border: 1px solid #c9a77c; color: #77532d; font-size: 10px; white-space: nowrap; }
.module-pagination { display: flex; align-items: center; justify-content: flex-end; gap: 8px; margin-top: 10px; }
.module-pagination button { padding: 5px 9px; }
.module-open-button { margin-left: 7px; padding: 4px 7px !important; font-size: 12px; }
.module-row { cursor: pointer; }
.module-row:hover { background: var(--hover); }
body.module-explorer-open { overflow: hidden; }
.artifact-stage-tabs {
  display: flex; flex-wrap: wrap; gap: 0; margin: 14px 0 9px;
  border-bottom: 1px solid var(--line-strong);
}
.artifact-stage-tabs button {
  border: 0; border-bottom: 3px solid transparent; background: transparent;
  padding: 9px 15px 8px; color: var(--text-secondary); font-weight: 600;
}
.artifact-stage-tabs button:hover { background: var(--hover); color: var(--text); }
.artifact-stage-tabs button.active {
  color: var(--text); border-bottom-color: var(--text); background: var(--selected);
}
.artifact-removal-tabs {
  display: flex; flex-wrap: wrap; gap: 6px; margin: 0 0 8px;
  padding: 7px 0; border-bottom: 1px solid var(--line-subtle);
}
.artifact-removal-tabs::before {
  content: 'Outlier removal'; align-self: center; margin-right: 4px;
  color: var(--text-muted); font-size: 12px; font-weight: 600;
}
.artifact-removal-tabs button {
  border: 1px solid var(--line); background: var(--surface); padding: 5px 10px;
  color: var(--text-secondary); font-size: 12px; font-weight: 600;
}
.artifact-removal-tabs button:hover { background: var(--hover); color: var(--text); }
.artifact-removal-tabs button.active {
  background: var(--selected); color: var(--text); border-color: var(--line-strong);
}
.artifact-tabs[data-stage-grouped="true"] {
  margin-top: 0; padding: 4px 0 9px; border-bottom: 1px solid var(--line-subtle);
}
.artifact-tabs[data-stage-grouped="true"] button { font-size: 13px; }
@media(max-width:760px) {
  .module-explorer { width: 100vw; padding: 18px 15px 32px; }
  .module-explorer-controls { grid-template-columns: 1fr; }
  .module-member-table { max-height: calc(100vh - 480px); }
  .artifact-stage-tabs button { flex: 1 1 50%; padding-inline: 8px; }
}
"""


MODULE_EXPLORER_HTML = r'''
<div class="module-explorer-overlay" id="moduleExplorerOverlay" hidden>
  <aside class="module-explorer" role="dialog" aria-modal="true" aria-labelledby="moduleExplorerTitle">
    <div class="module-explorer-head">
      <div><h2 id="moduleExplorerTitle">Module Explorer</h2><div class="sub" id="moduleExplorerRecord"></div></div>
      <button type="button" id="moduleExplorerClose">Close</button>
    </div>
    <div class="module-explorer-controls">
      <select id="moduleExplorerModule" aria-label="Module"></select>
      <input id="moduleExplorerSearch" type="search" placeholder="Search feature, gene, PROBEID, UniProt">
      <select id="moduleExplorerSort" aria-label="Sort members">
        <option value="rank">Hub rank (ascending)</option>
        <option value="abs_desc">|kME| (descending)</option>
        <option value="signed_desc">Signed kME (descending)</option>
        <option value="feature">Feature ID (A-Z)</option>
      </select>
      <button type="button" id="moduleExplorerExport">Export CSV</button>
    </div>
    <div class="module-explorer-summary" id="moduleExplorerSummary"></div>
    <div class="module-explorer-note">kME is within-network module-membership evidence, not a p-value or disease association. Raw-unassigned features remain visible but their undefined kME and rank are shown as NA.</div>
    <div class="table-wrap module-member-table">
      <table><thead><tr><th>Hub rank</th><th>Feature ID</th><th>Gene ID</th><th>Symbol</th><th>PROBEID</th><th>UniProt</th><th>Signed kME</th><th>|kME|</th><th>Status</th></tr></thead><tbody id="moduleExplorerRows"></tbody></table>
    </div>
    <div class="module-pagination"><button type="button" id="moduleExplorerPrev">Previous</button><span id="moduleExplorerPage"></span><button type="button" id="moduleExplorerNext">Next</button></div>
  </aside>
</div>
'''


MODULE_EXPLORER_JS = r'''
const moduleMemberRecords=(window.MODULE_MEMBER_INDEX&&window.MODULE_MEMBER_INDEX.records)||{};
const moduleExplorerState={recordId:null,moduleColor:null,payload:null,query:'',sort:'rank',page:1,pageSize:100,loadingToken:0};
const moduleShardPromises=new Map();
function moduleEntry(recordId){return moduleMemberRecords[recordId]||null}
function moduleKey(recordId,moduleColor){return `${recordId}::${moduleColor}`}
function moduleDefaultColor(entry){
  const modules=Object.entries(entry?.modules||{}).sort((a,b)=>Number(a[1].fully_raw_unassigned)-Number(b[1].fully_raw_unassigned)||b[1].count-a[1].count||a[0].localeCompare(b[0]));
  return modules[0]?.[0]||null;
}
function loadModuleShard(recordId,moduleColor){
  const entry=moduleEntry(recordId),meta=entry?.modules?.[moduleColor],key=moduleKey(recordId,moduleColor);
  if(!meta)return Promise.reject(new Error('Module is not indexed'));
  if(window.MODULE_MEMBER_SHARDS?.[key])return Promise.resolve(window.MODULE_MEMBER_SHARDS[key]);
  if(moduleShardPromises.has(key))return moduleShardPromises.get(key);
  const promise=new Promise((resolve,reject)=>{const script=document.createElement('script');script.src=meta.path;script.onload=()=>window.MODULE_MEMBER_SHARDS?.[key]?resolve(window.MODULE_MEMBER_SHARDS[key]):reject(new Error('Module shard did not register'));script.onerror=()=>reject(new Error(`Could not load ${meta.path}`));document.head.appendChild(script)});
  moduleShardPromises.set(key,promise);return promise;
}
function moduleNumber(v,d=3){return v===null||v===undefined||v===''?'NA':Number(v).toLocaleString(undefined,{maximumFractionDigits:d})}
function moduleFilteredRows(){
  if(!moduleExplorerState.payload)return [];
  const q=moduleExplorerState.query.trim().toLowerCase();let rows=moduleExplorerState.payload.rows.filter(r=>!q||r.slice(0,5).some(v=>String(v||'').toLowerCase().includes(q)));
  rows=[...rows];const nullLast=(v,fallback)=>v===null||v===undefined?fallback:Number(v);
  if(moduleExplorerState.sort==='rank')rows.sort((a,b)=>nullLast(a[7],Infinity)-nullLast(b[7],Infinity)||nullLast(b[6],-Infinity)-nullLast(a[6],-Infinity)||String(a[0]).localeCompare(String(b[0])));
  else if(moduleExplorerState.sort==='abs_desc')rows.sort((a,b)=>nullLast(b[6],-Infinity)-nullLast(a[6],-Infinity)||String(a[0]).localeCompare(String(b[0])));
  else if(moduleExplorerState.sort==='signed_desc')rows.sort((a,b)=>nullLast(b[5],-Infinity)-nullLast(a[5],-Infinity)||String(a[0]).localeCompare(String(b[0])));
  else rows.sort((a,b)=>String(a[0]).localeCompare(String(b[0])));
  return rows;
}
function renderModuleExplorer(){
  const entry=moduleEntry(moduleExplorerState.recordId),meta=entry?.modules?.[moduleExplorerState.moduleColor];if(!entry||!meta)return;
  document.getElementById('moduleExplorerRecord').textContent=`${entry.title||moduleExplorerState.recordId} · ${entry.network_id||''}`;
  const select=document.getElementById('moduleExplorerModule');select.innerHTML=Object.entries(entry.modules).sort((a,b)=>Number(a[1].fully_raw_unassigned)-Number(b[1].fully_raw_unassigned)||b[1].count-a[1].count||a[0].localeCompare(b[0])).map(([color,m])=>`<option value="${esc(color)}" ${color===moduleExplorerState.moduleColor?'selected':''}>${esc(color)} · ${m.count.toLocaleString()} features${m.fully_raw_unassigned?' · raw-unassigned':''}</option>`).join('');
  const rows=moduleFilteredRows(),pages=Math.max(1,Math.ceil(rows.length/moduleExplorerState.pageSize));moduleExplorerState.page=Math.min(Math.max(1,moduleExplorerState.page),pages);const start=(moduleExplorerState.page-1)*moduleExplorerState.pageSize,shown=rows.slice(start,start+moduleExplorerState.pageSize);
  document.getElementById('moduleExplorerSummary').innerHTML=`<span><span class="swatch" style="background:${esc(moduleExplorerState.moduleColor)}"></span><b>${esc(moduleExplorerState.moduleColor)}</b></span><span><b>${meta.count.toLocaleString()}</b> total features</span><span><b>${meta.assigned_kme_count.toLocaleString()}</b> with kME</span><span><b>${meta.raw_unassigned_count.toLocaleString()}</b> raw-unassigned</span><span>median |kME| <b>${moduleNumber(meta.median_abs_kme)}</b></span><span>max |kME| <b>${moduleNumber(meta.max_abs_kme)}</b></span><span>${rows.length.toLocaleString()} search matches</span>`;
  document.getElementById('moduleExplorerRows').innerHTML=shown.map(r=>`<tr class="${r[8]?'raw-unassigned':''}"><td class="num">${moduleNumber(r[7],0)}</td><td>${esc(r[0])}</td><td>${esc(r[1])}</td><td>${esc(r[2])}</td><td>${esc(r[3])}</td><td>${esc(r[4])}</td><td class="num">${moduleNumber(r[5])}</td><td class="num">${moduleNumber(r[6])}</td><td>${r[8]?'<span class="raw-badge">raw-unassigned</span>':'assigned'}</td></tr>`).join('')||'<tr><td colspan="9" class="empty">No matching features</td></tr>';
  document.getElementById('moduleExplorerPage').textContent=`${start+1}-${Math.min(start+shown.length,rows.length)} of ${rows.length.toLocaleString()} · page ${moduleExplorerState.page}/${pages}`;
  document.getElementById('moduleExplorerPrev').disabled=moduleExplorerState.page<=1;document.getElementById('moduleExplorerNext').disabled=moduleExplorerState.page>=pages;
}
async function selectModuleMembers(recordId,moduleColor){
  const token=++moduleExplorerState.loadingToken;moduleExplorerState.recordId=recordId;moduleExplorerState.moduleColor=moduleColor;moduleExplorerState.payload=null;moduleExplorerState.page=1;
  document.getElementById('moduleExplorerRows').innerHTML='<tr><td colspan="9" class="empty">Loading module members…</td></tr>';
  try{const payload=await loadModuleShard(recordId,moduleColor);if(token!==moduleExplorerState.loadingToken)return;moduleExplorerState.payload=payload;renderModuleExplorer()}catch(error){if(token!==moduleExplorerState.loadingToken)return;document.getElementById('moduleExplorerRows').innerHTML=`<tr><td colspan="9" class="empty">${esc(error.message)}</td></tr>`}
}
function openModuleExplorer(recordId,moduleColor){
  const entry=moduleEntry(recordId);if(!entry)return;const color=entry.modules[moduleColor]?moduleColor:moduleDefaultColor(entry);if(!color)return;
  document.getElementById('moduleExplorerOverlay').hidden=false;document.body.classList.add('module-explorer-open');moduleExplorerState.query='';moduleExplorerState.sort='rank';document.getElementById('moduleExplorerSearch').value='';document.getElementById('moduleExplorerSort').value='rank';selectModuleMembers(recordId,color);
}
function closeModuleExplorer(){document.getElementById('moduleExplorerOverlay').hidden=true;document.body.classList.remove('module-explorer-open')}
window.openModuleExplorer=openModuleExplorer;
function moduleCsvCell(v){const s=String(v??'');return /[",\n]/.test(s)?`"${s.replaceAll('"','""')}"`:s}
function exportModuleMembers(){
  const rows=moduleFilteredRows(),head=['hub_rank','feature_id','gene_id','symbol','PROBEID','UNIPROT','signed_kME','abs_kME','is_raw_unassigned'];const data=rows.map(r=>[r[7],r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[8]]);const blob=new Blob(['\ufeff'+[head,...data].map(row=>row.map(moduleCsvCell).join(',')).join('\n')],{type:'text/csv;charset=utf-8'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=`${moduleExplorerState.recordId.replace(/[^a-zA-Z0-9_-]+/g,'_')}__${moduleExplorerState.moduleColor}__members.csv`;a.click();URL.revokeObjectURL(a.href);
}
document.getElementById('moduleExplorerClose').onclick=closeModuleExplorer;document.getElementById('moduleExplorerOverlay').onclick=e=>{if(e.target.id==='moduleExplorerOverlay')closeModuleExplorer()};document.addEventListener('keydown',e=>{if(e.key==='Escape'&&!document.getElementById('moduleExplorerOverlay').hidden)closeModuleExplorer()});
document.getElementById('moduleExplorerModule').onchange=e=>selectModuleMembers(moduleExplorerState.recordId,e.target.value);document.getElementById('moduleExplorerSearch').oninput=e=>{moduleExplorerState.query=e.target.value;moduleExplorerState.page=1;renderModuleExplorer()};document.getElementById('moduleExplorerSort').onchange=e=>{moduleExplorerState.sort=e.target.value;moduleExplorerState.page=1;renderModuleExplorer()};document.getElementById('moduleExplorerPrev').onclick=()=>{moduleExplorerState.page--;renderModuleExplorer()};document.getElementById('moduleExplorerNext').onclick=()=>{moduleExplorerState.page++;renderModuleExplorer()};document.getElementById('moduleExplorerExport').onclick=exportModuleMembers;
document.addEventListener('click',e=>{const moduleButton=e.target.closest('[data-open-modules]');if(moduleButton){openModuleExplorer(moduleButton.dataset.openModules,moduleButton.dataset.moduleColor||null);return}const row=e.target.closest('[data-module-record][data-module-color]');if(row)openModuleExplorer(row.dataset.moduleRecord,row.dataset.moduleColor)});
const artifactStageDefinitions=[
  {id:'before',label:'Before Preprocess',artifacts:['feature_selection','before_heatmap','before_pca','before_pca_3d','after_outlier_heatmap','after_outlier_pca','after_outlier_pca_3d']},
  {id:'preprocess',label:'Preprocess',artifacts:['cut_diagnostic','selected_cut','after_outlier_cut_diagnostic']},
  {id:'after',label:'After Preprocess',artifacts:['after_heatmap','after_pca','after_pca_3d','after_outlier_after_heatmap','after_outlier_after_pca','after_outlier_after_pca_3d']},
  {id:'modules',label:'Find Modules',artifacts:['soft_threshold','gene_dendrogram','me_dendrogram','module_sizes','after_outlier_soft_threshold','after_outlier_gene_dendrogram','after_outlier_me_dendrogram','after_outlier_module_sizes']},
  {id:'downstream',label:'Downstream',artifacts:['diagnosis_pairwise_hd','diagnosis_heatmap','cross_omics_similarity','after_outlier_diagnosis_pairwise_hd','after_outlier_diagnosis_heatmap','after_outlier_cross_omics_similarity']},
  {id:'other',label:'Other',artifacts:[]},
];
const artifactStageById=new Map(artifactStageDefinitions.flatMap(stage=>stage.artifacts.map(id=>[id,stage.id])));
const artifactStageState={};
const artifactRemovalState={};
const artifactRemovalDefinitions=[
  {id:'before',label:'Before Outlier Removal',artifacts:['feature_selection','before_heatmap','before_pca','before_pca_3d','cut_diagnostic','selected_cut','after_heatmap','after_pca','after_pca_3d','soft_threshold','gene_dendrogram','me_dendrogram','module_sizes','diagnosis_pairwise_hd','diagnosis_heatmap','cross_omics_similarity']},
  {id:'after',label:'After Outlier Removal',artifacts:['after_outlier_heatmap','after_outlier_pca','after_outlier_pca_3d','after_outlier_cut_diagnostic','after_outlier_after_heatmap','after_outlier_after_pca','after_outlier_after_pca_3d','after_outlier_soft_threshold','after_outlier_gene_dendrogram','after_outlier_me_dendrogram','after_outlier_module_sizes','after_outlier_diagnosis_pairwise_hd','after_outlier_diagnosis_heatmap','after_outlier_cross_omics_similarity']},
];
const artifactRemovalById=new Map(artifactRemovalDefinitions.flatMap(group=>group.artifacts.map(id=>[id,group.id])));
function artifactStageOf(id){return artifactStageById.get(id)||'other'}
function artifactRemovalOf(id){return artifactRemovalById.get(id)||null}
function artifactRemovalPairKey(id){return String(id||'').replace(/^after_outlier_/,'').replace(/^before_/,'')}
function decorateArtifactStageBar(tabBar,stateKey,buttonSelector){
  if(!tabBar||tabBar.dataset.stageGrouped==='true')return;
  const buttons=[...tabBar.querySelectorAll(buttonSelector)];if(!buttons.length)return;
  const available=artifactStageDefinitions.filter(stage=>buttons.some(button=>artifactStageOf(button.dataset.artifact||button.dataset.comparisonArtifact)===stage.id));
  if(!available.length)return;
  const active=buttons.find(button=>button.classList.contains('active')),activeStage=artifactStageOf(active?.dataset.artifact||active?.dataset.comparisonArtifact);
  let chosen=artifactStageState[stateKey];if(!available.some(stage=>stage.id===chosen))chosen=available.some(stage=>stage.id===activeStage)?activeStage:available[0].id;artifactStageState[stateKey]=chosen;
  const nav=document.createElement('div');nav.className='artifact-stage-tabs';nav.setAttribute('aria-label','Analysis stage');
  for(const stage of available){const button=document.createElement('button');button.type='button';button.dataset.artifactStage=stage.id;button.className=stage.id===chosen?'active':'';button.textContent=stage.label;button.onclick=()=>{artifactStageState[stateKey]=stage.id;const first=buttons.find(item=>artifactStageOf(item.dataset.artifact||item.dataset.comparisonArtifact)===stage.id);if(first)first.click()};nav.appendChild(button)}
  const stageButtons=buttons.filter(button=>artifactStageOf(button.dataset.artifact||button.dataset.comparisonArtifact)===chosen);
  const removalGroups=artifactRemovalDefinitions.filter(group=>stageButtons.some(button=>artifactRemovalOf(button.dataset.artifact||button.dataset.comparisonArtifact)===group.id));
  let removalChosen=null,removalNav=null;
  if(removalGroups.length>1){const activeId=active?.dataset.artifact||active?.dataset.comparisonArtifact;const activeRemoval=artifactRemovalOf(activeId);removalChosen=artifactRemovalState[stateKey];if(!removalGroups.some(group=>group.id===removalChosen))removalChosen=removalGroups.some(group=>group.id===activeRemoval)?activeRemoval:removalGroups[0].id;artifactRemovalState[stateKey]=removalChosen;removalNav=document.createElement('div');removalNav.className='artifact-removal-tabs';removalNav.setAttribute('aria-label','Outlier removal state');for(const group of removalGroups){const button=document.createElement('button');button.type='button';button.dataset.artifactRemoval=group.id;button.className=group.id===removalChosen?'active':'';button.textContent=group.label;button.onclick=()=>{artifactRemovalState[stateKey]=group.id;const pairKey=artifactRemovalPairKey(activeId);const counterpart=stageButtons.find(item=>{const id=item.dataset.artifact||item.dataset.comparisonArtifact;return artifactRemovalOf(id)===group.id&&artifactRemovalPairKey(id)===pairKey});const fallback=stageButtons.find(item=>artifactRemovalOf(item.dataset.artifact||item.dataset.comparisonArtifact)===group.id);(counterpart||fallback)?.click()};removalNav.appendChild(button)}}
  for(const button of buttons){const id=button.dataset.artifact||button.dataset.comparisonArtifact;const stageVisible=artifactStageOf(id)===chosen;const removalVisible=!removalChosen||artifactRemovalOf(id)===removalChosen;button.style.display=stageVisible&&removalVisible?'':'none'}
  tabBar.dataset.stageGrouped='true';tabBar.parentNode.insertBefore(nav,tabBar);if(removalNav)tabBar.parentNode.insertBefore(removalNav,tabBar);
}
function decorateArtifactStages(){
  const comparisonBar=document.querySelector('#cellComparison > .artifact-tabs');if(comparisonBar){const key=`comparison:${document.getElementById('source')?.value||''}:${document.getElementById('cell')?.value||''}`;decorateArtifactStageBar(comparisonBar,key,'[data-comparison-artifact]')}
  const detailBar=document.querySelector('#detail > .artifact-tabs');if(detailBar)decorateArtifactStageBar(detailBar,`detail:${selectedId}`,'[data-artifact]');
}
let moduleDecoratePending=false;
function decorateModuleEntrypoints(){
  moduleDecoratePending=false;
  decorateArtifactStages();
  document.querySelectorAll('.condition-panel').forEach(panel=>{const inspect=panel.querySelector('[data-open-record]'),id=inspect?.dataset.openRecord;if(!id||!moduleEntry(id)||panel.querySelector('[data-open-modules]'))return;const button=document.createElement('button');button.type='button';button.className='module-open-button';button.dataset.openModules=id;button.textContent='Modules & kME';panel.querySelector('.condition-links')?.appendChild(button)});
  const detailActions=document.querySelector('#detail .detail-head .actions');if(detailActions&&moduleEntry(selectedId)&&!detailActions.querySelector('[data-open-modules]')){const button=document.createElement('button');button.type='button';button.dataset.openModules=selectedId;button.textContent='Modules & kME';detailActions.appendChild(button)}
  document.querySelectorAll('#detail .lower > div:first-child tbody tr').forEach(row=>{if(row.dataset.moduleRecord)return;const cell=row.querySelector('td'),color=cell?.textContent.trim();if(!color||!moduleEntry(selectedId)?.modules?.[color])return;row.classList.add('module-row');row.dataset.moduleRecord=selectedId;row.dataset.moduleColor=color;row.title='Open module members and kME ranks';const button=document.createElement('button');button.type='button';button.className='module-open-button';button.dataset.openModules=selectedId;button.dataset.moduleColor=color;button.textContent='View members';cell.appendChild(button)});
  document.querySelectorAll('#compareGrid .compare-card').forEach((card,i)=>{const id=compareIds[i];if(!id||!moduleEntry(id)||card.querySelector('[data-open-modules]'))return;const button=document.createElement('button');button.type='button';button.className='module-open-button';button.dataset.openModules=id;button.textContent='Modules & kME';card.appendChild(button)});
}
function scheduleModuleDecoration(){if(moduleDecoratePending)return;moduleDecoratePending=true;queueMicrotask(decorateModuleEntrypoints)}
new MutationObserver(scheduleModuleDecoration).observe(document.getElementById('cellComparison').parentElement,{childList:true,subtree:true});scheduleModuleDecoration();
const moduleParams=new URLSearchParams(location.search);const deepRecord=moduleParams.get('module_record');if(deepRecord&&moduleEntry(deepRecord))setTimeout(()=>openModuleExplorer(deepRecord,moduleParams.get('module_color')),0);
'''


CROSS_OMICS_EXPLORER_CSS = r"""
.cross-omics-overlay[hidden] { display: none; }
.cross-omics-overlay {
  position: fixed; inset: 0; z-index: 1010; background: rgba(20,20,20,.38);
  display: flex; justify-content: flex-end;
}
.cross-omics-explorer {
  width: min(1080px, 97vw); height: 100vh; overflow: auto;
  background: var(--surface); border-left: 1px solid var(--line-strong);
  box-shadow: -12px 0 36px rgba(0,0,0,.12); padding: 24px 28px 40px;
}
.cross-omics-head { display: flex; gap: 20px; align-items: flex-start; border-bottom: 1px solid var(--line); padding-bottom: 15px; }
.cross-omics-head h2 { margin: 0 0 3px; font-size: 23px; }
.cross-omics-head button { margin-left: auto; padding: 7px 12px; }
.cross-omics-controls { display: grid; grid-template-columns: minmax(170px,.7fr) minmax(170px,.7fr) minmax(210px,1fr) minmax(160px,.7fr) auto; gap: 8px; margin: 16px 0 10px; }
.cross-omics-controls select, .cross-omics-controls input, .cross-omics-controls button { min-width: 0; padding: 8px 9px; }
.cross-omics-heatmap-section { margin: 17px 0 13px; }
.cross-omics-heatmap-section h3 { margin: 0 0 4px; font-size: 16px; }
.cross-omics-heatmap-section .sub { margin: 0 0 9px; }
.cross-omics-heatmap { max-height: 36vh; overflow: auto; }
.cross-omics-heatmap table { border-collapse: separate; border-spacing: 0; }
.cross-omics-heatmap th, .cross-omics-heatmap td { text-align: center; white-space: nowrap; }
.cross-omics-heatmap thead th { position: sticky; top: 0; z-index: 2; background: var(--surface); }
.cross-omics-heatmap th[scope=row] { position: sticky; left: 0; z-index: 1; background: var(--surface); }
.cross-omics-heatmap thead th:first-child { left: 0; z-index: 3; }
.cross-omics-heatmap td { padding: 2px; }
.cross-heatmap-cell { width: 100%; min-width: 4rem; border: 1px solid transparent; padding: 5px; font-variant-numeric: tabular-nums; }
.cross-heatmap-cell[aria-pressed=true] { outline: 2px solid #b85017; outline-offset: -2px; font-weight: bold; }
.cross-heatmap-cell:focus-visible { outline: 2px solid #101820; outline-offset: -4px; }
.cross-omics-summary { display: flex; flex-wrap: wrap; gap: 7px 16px; padding: 10px 12px; background: var(--panel); border: 1px solid var(--line); font-size: 13px; }
.cross-gene-groups { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 8px; margin: 10px 0; }
.cross-gene-groups details { border: 1px solid var(--line); background: var(--surface); padding: 8px 10px; min-width: 0; }
.cross-gene-groups summary { cursor: pointer; font-weight: 600; }
.cross-gene-list { max-height: 145px; overflow: auto; overflow-wrap: anywhere; margin: 8px 0 0; line-height: 1.45; color: var(--text-muted); }
.forest-track { position: relative; min-width: 170px; height: 27px; background: var(--panel); }
.forest-track::before { content: ''; position: absolute; inset: 0 auto 0 50%; border-left: 1px dashed #6b7280; }
.forest-interval, .forest-point { position: absolute; top: 50%; background: var(--forest-color); }
.forest-interval { height: 2px; transform: translateY(-50%); }
.forest-point { width: 9px; height: 9px; transform: translate(-50%, -50%); }
.forest-rna { --forest-color: #216683; }
.forest-protein { --forest-color: #87532d; }
.forest-axis { display: flex; justify-content: space-between; font-variant-numeric: tabular-nums; }
.cross-omics-diagnosis { margin-top: 22px; }
.cross-omics-diagnosis .overlap-forest { max-height: 430px; }
.cross-omics-note { color: var(--text-muted); font-size: 12px; margin: 9px 0 13px; }
.cross-omics-table { max-height: calc(100vh - 350px); }
.cross-omics-table td.num { text-align: right; font-variant-numeric: tabular-nums; }
.cross-best-badge { display: inline-block; padding: 1px 5px; border: 1px solid var(--line-strong); font-size: 10px; }
.cross-overlap-button { margin-left: 7px; padding: 4px 7px !important; font-size: 12px; }
body.cross-omics-open { overflow: hidden; }
@media(max-width:820px) {
  .cross-omics-explorer { width: 100vw; padding: 18px 15px 32px; }
  .cross-omics-controls { grid-template-columns: 1fr; }
  .cross-gene-groups { grid-template-columns: 1fr; }
  .cross-omics-table { max-height: calc(100vh - 540px); }
}
"""


CROSS_OMICS_EXPLORER_HTML = r'''
<div class="cross-omics-overlay" id="crossOmicsOverlay" hidden>
  <aside class="cross-omics-explorer" role="dialog" aria-modal="true" aria-labelledby="crossOmicsTitle">
    <div class="cross-omics-head">
      <div><h2 id="crossOmicsTitle">RNA–protein module comparison</h2><div class="sub" id="crossOmicsRecord"></div></div>
      <button type="button" id="crossOmicsClose">Close</button>
    </div>
    <section class="cross-omics-heatmap-section" aria-labelledby="crossOmicsHeatmapTitle">
      <h3 id="crossOmicsHeatmapTitle">kME-weighted Jaccard similarity</h3>
      <p class="sub">Select an RNA row × Protein column cell. Darker cells indicate higher similarity. • RNA best match · ★ reciprocal best</p>
      <div class="table-wrap cross-omics-heatmap" id="crossOmicsHeatmap"></div>
    </section>
    <div class="cross-omics-controls">
      <select id="crossTranscriptModule" aria-label="Transcriptomics module"></select>
      <select id="crossProteomicsModule" aria-label="Proteomics module"></select>
      <input id="crossOmicsSearch" type="search" placeholder="Search overlapping gene">
      <select id="crossOmicsSort" aria-label="Sort overlapping genes">
        <option value="intersection_desc">Shared kME contribution (descending)</option>
        <option value="transcript_desc">Transcriptomics |kME| (descending)</option>
        <option value="proteomics_desc">Proteomics |kME| (descending)</option>
        <option value="gene">Gene (A-Z)</option>
      </select>
      <button type="button" id="crossOmicsExport">Export CSV</button>
    </div>
    <div class="cross-omics-summary" id="crossOmicsSummary"></div>
    <h3>Module genes</h3>
    <p class="sub">Lists use mapped genes in each module. The heatmap scores use the shared measurement universe.</p>
    <div class="cross-gene-groups" id="crossOmicsGeneGroups" aria-live="polite"></div>
    <section class="cross-omics-diagnosis" aria-labelledby="crossOmicsDiagnosisTitle">
      <h3 id="crossOmicsDiagnosisTitle">Module-diagnosis · disease vs HD</h3>
      <div id="crossOmicsDiagnosis" aria-live="polite"></div>
    </section>
    <h3>Shared gene kME detail</h3>
    <div class="cross-omics-note">Only genes present in both selected modules are listed below. Weighted Jaccard also includes genes found in only one module in its denominator.</div>
    <div class="table-wrap cross-omics-table">
      <table><thead><tr><th>Gene</th><th>Transcriptomics |kME|</th><th>Proteomics |kME|</th><th>min(|kME|)</th><th>max(|kME|)</th></tr></thead><tbody id="crossOmicsRows"></tbody></table>
    </div>
    <div class="module-pagination"><button type="button" id="crossOmicsPrev">Previous</button><span id="crossOmicsPage"></span><button type="button" id="crossOmicsNext">Next</button></div>
  </aside>
</div>
'''


CROSS_OMICS_EXPLORER_JS = r'''
function crossRecordIsOtherCell(recordId,entry){const base=String(recordId).replace(/::after_outlier$/,'');const parts=base.split('::');return String(entry?.cell_type||'').trim().toLowerCase()==='other'||String(parts[parts.length-1]||'').trim().toLowerCase()==='other'}
const crossOmicsView=window.CROSS_OMICS_VIEW||{diagnosis:{}};
const crossOverlapRecords=Object.fromEntries(Object.entries((window.CROSS_OMICS_OVERLAP_INDEX&&window.CROSS_OMICS_OVERLAP_INDEX.records)||{}).filter(([recordId,entry])=>!crossRecordIsOtherCell(recordId,entry)));
const crossOverlapState={recordId:null,transcriptModule:null,proteomicsModule:null,payload:null,transcriptMembers:null,proteomicsMembers:null,geneGroupError:null,query:'',sort:'intersection_desc',page:1,pageSize:100,loadingToken:0};
const crossShardPromises=new Map();
function crossEntry(recordId){return crossOverlapRecords[recordId]||null}
function crossVariantRecordId(recordId,artifactId){return String(artifactId||'').startsWith('after_outlier_')?`${recordId}::after_outlier`:recordId}
function crossPairObject(entry,row){return Object.fromEntries(entry.pair_columns.map((column,index)=>[column,row[index]]))}
function crossPairs(entry){return (entry?.pairs||[]).map(row=>crossPairObject(entry,row))}
function crossSelectedPair(){const entry=crossEntry(crossOverlapState.recordId);return crossPairs(entry).find(pair=>pair.transcript_module===crossOverlapState.transcriptModule&&pair.proteomics_module===crossOverlapState.proteomicsModule)||null}
function loadCrossShard(recordId){
  const entry=crossEntry(recordId);if(!entry)return Promise.reject(new Error('Cross-omics overlap is not indexed'));
  if(window.CROSS_OMICS_OVERLAP_SHARDS?.[recordId])return Promise.resolve(window.CROSS_OMICS_OVERLAP_SHARDS[recordId]);
  if(crossShardPromises.has(recordId))return crossShardPromises.get(recordId);
  const promise=new Promise((resolve,reject)=>{const script=document.createElement('script');script.src=entry.path;script.onload=()=>window.CROSS_OMICS_OVERLAP_SHARDS?.[recordId]?resolve(window.CROSS_OMICS_OVERLAP_SHARDS[recordId]):reject(new Error('Overlap shard did not register'));script.onerror=()=>reject(new Error(`Could not load ${entry.path}`));document.head.appendChild(script)});crossShardPromises.set(recordId,promise);return promise;
}
function crossNumber(value,digits=3){return value===null||value===undefined||value===''?'NA':Number(value).toLocaleString(undefined,{maximumFractionDigits:digits})}
function crossFilteredRows(){
  if(!crossOverlapState.payload)return [];const q=crossOverlapState.query.trim().toLowerCase();let rows=crossOverlapState.payload.rows.filter(row=>row[0]===crossOverlapState.transcriptModule&&row[1]===crossOverlapState.proteomicsModule&&(!q||String(row[2]).toLowerCase().includes(q)));rows=[...rows];
  const number=(v,fallback)=>v===null||v===undefined?fallback:Number(v);
  if(crossOverlapState.sort==='intersection_desc')rows.sort((a,b)=>number(b[5],-Infinity)-number(a[5],-Infinity)||String(a[2]).localeCompare(String(b[2])));
  else if(crossOverlapState.sort==='transcript_desc')rows.sort((a,b)=>number(b[3],-Infinity)-number(a[3],-Infinity)||String(a[2]).localeCompare(String(b[2])));
  else if(crossOverlapState.sort==='proteomics_desc')rows.sort((a,b)=>number(b[4],-Infinity)-number(a[4],-Infinity)||String(a[2]).localeCompare(String(b[2])));
  else rows.sort((a,b)=>String(a[2]).localeCompare(String(b[2])));return rows;
}
function updateCrossModuleOptions(){
  const entry=crossEntry(crossOverlapState.recordId);if(!entry)return;
  const transcript=document.getElementById('crossTranscriptModule');transcript.innerHTML=entry.transcript_modules.map(module=>`<option value="${esc(module)}" ${module===crossOverlapState.transcriptModule?'selected':''}>RNA · ${esc(module)}</option>`).join('');
  const proteomics=document.getElementById('crossProteomicsModule');proteomics.innerHTML=entry.proteomics_modules.map(module=>`<option value="${esc(module)}" ${module===crossOverlapState.proteomicsModule?'selected':''}>Protein · ${esc(module)}</option>`).join('');
}
function renderCrossOmicsHeatmap(){
  const entry=crossEntry(crossOverlapState.recordId),container=document.getElementById('crossOmicsHeatmap');if(!entry||!container)return;
  const pairs=crossPairs(entry),pairByModules=new Map(pairs.map(pair=>[`${pair.transcript_module}\u0000${pair.proteomics_module}`,pair])),maxWeighted=Math.max(0,...pairs.map(pair=>Number(pair.weighted_jaccard)||0)),scrollTop=container.scrollTop,scrollLeft=container.scrollLeft;
  const header=entry.proteomics_modules.map(module=>`<th scope="col">${esc(module)}</th>`).join('');
  const rows=entry.transcript_modules.map(transcriptModule=>`<tr><th scope="row">${esc(transcriptModule)}</th>${entry.proteomics_modules.map(proteomicsModule=>{
    const pair=pairByModules.get(`${transcriptModule}\u0000${proteomicsModule}`);if(!pair)return '<td>NA</td>';
    const weighted=Number(pair.weighted_jaccard)||0,intensity=maxWeighted>0?weighted/maxWeighted:0,alpha=(.05+.9*intensity).toFixed(3),marker=pair.reciprocal_best?' ★':pair.transcript_best?' •':'',selected=transcriptModule===crossOverlapState.transcriptModule&&proteomicsModule===crossOverlapState.proteomicsModule;
    return `<td><button type="button" class="cross-heatmap-cell" data-cross-rna="${esc(transcriptModule)}" data-cross-protein="${esc(proteomicsModule)}" aria-pressed="${selected}" aria-label="RNA ${esc(transcriptModule)}, Protein ${esc(proteomicsModule)}, weighted Jaccard ${weighted.toFixed(4)}${pair.reciprocal_best?', reciprocal best':pair.transcript_best?', RNA best match':''}" style="background:rgba(31,100,133,${alpha});color:${intensity>.72?'white':'#16232d'}">${weighted.toFixed(3)}${marker}</button></td>`;
  }).join('')}</tr>`).join('');
  container.innerHTML=`<table aria-label="kME-weighted Jaccard heatmap"><thead><tr><th scope="col">RNA ↓ / Protein →</th>${header}</tr></thead><tbody>${rows}</tbody></table>`;container.scrollTop=scrollTop;container.scrollLeft=scrollLeft;
}
function crossMemberGeneSet(payload,omic){
  const genes=new Set();if(!payload?.rows)return genes;
  for(const row of payload.rows){if(row[8])continue;const value=omic==='protein'?row[2]:(row[1]||row[2]||row[0]);const gene=String(value||'').trim().toUpperCase();if(!gene||gene==='---'||gene==='NA'||gene==='NAN'||/[|;,/]/.test(gene))continue;genes.add(gene)}return genes;
}
function renderCrossGeneGroups(){
  const container=document.getElementById('crossOmicsGeneGroups');if(!container)return;
  if(crossOverlapState.geneGroupError){container.innerHTML=`<div class="empty" style="grid-column:1/-1">${esc(crossOverlapState.geneGroupError)}</div>`;return}
  if(!crossOverlapState.transcriptMembers||!crossOverlapState.proteomicsMembers){container.innerHTML='<div class="empty" style="grid-column:1/-1">Loading shared, RNA-only, and Protein-only genes…</div>';return}
  const rna=crossMemberGeneSet(crossOverlapState.transcriptMembers,'rna'),protein=new Set(crossOverlapState.proteomicsMembers.map(gene=>String(gene).trim().toUpperCase()).filter(Boolean)),shared=[...rna].filter(gene=>protein.has(gene)).sort(),rnaOnly=[...rna].filter(gene=>!protein.has(gene)).sort(),proteinOnly=[...protein].filter(gene=>!rna.has(gene)).sort();
  const groups=[['Shared genes',shared],['RNA-only genes',rnaOnly],['Protein-only genes',proteinOnly]];
  container.innerHTML=groups.map(([label,genes])=>`<details open><summary>${label} (${genes.length.toLocaleString()})</summary><p class="cross-gene-list">${genes.length?genes.map(esc).join(', '):'None'}</p></details>`).join('');
}
function crossDiagnosisInterval(r,n){
  if(!Number.isFinite(r)||Math.abs(r)>=1||!Number.isFinite(n)||n<=3)return null;
  const z=Math.atanh(r),margin=1.959964/Math.sqrt(n-3);return [Math.tanh(z-margin),Math.tanh(z+margin)];
}
function renderCrossOmicsDiagnosis(){
  const container=document.getElementById('crossOmicsDiagnosis'),entry=crossEntry(crossOverlapState.recordId),pair=crossSelectedPair();if(!container||!entry||!pair)return;
  const proteinRecordId=`proteomics::${entry.specimen}::${entry.proteomics_condition}`,diagnosis=crossOmicsView.diagnosis||{},modules=[['RNA',pair.transcript_module,diagnosis[crossOverlapState.recordId]?.[pair.transcript_module]||[]],['Protein',pair.proteomics_module,diagnosis[proteinRecordId]?.[pair.proteomics_module]||[]]],contrasts=[...new Set(modules.flatMap(([, ,rows])=>rows.map(row=>row[0])))];
  if(!contrasts.length){container.textContent='No saved diagnosis statistics for these modules.';return}
  const rows=contrasts.flatMap(contrast=>modules.map(([assay,color,results])=>{const result=results.find(row=>row[0]===contrast),[,r,p,q,nHD,nTarget,n,ready,reason]=result||[];return {contrast,assay,color,r,p,q,nHD,nTarget,ready,reason,interval:crossDiagnosisInterval(r,n)}})),farthest=Math.max(.25,...rows.flatMap(({r,interval})=>[Number.isFinite(r)?Math.abs(r):0,...(interval||[]).map(Math.abs)])),limit=Math.ceil(farthest*4)/4,position=value=>50*(1+value/limit);
  const body=rows.map(({contrast,assay,color,r,p,q,nHD,nTarget,ready,reason,interval})=>{const range=interval?`[${interval.map(value=>value.toFixed(3)).join(', ')}]`:'—',line=interval?`<span class="forest-interval" style="left:${position(interval[0])}%;width:${position(interval[1])-position(interval[0])}%"></span>`:'',point=Number.isFinite(r)?`<span class="forest-point" style="left:${position(r)}%"></span>`:'',status=ready?'':reason||(Number.isFinite(r)?'Low sample support':'Unavailable');return `<tr><th scope="row">${esc(contrast)}</th><td>${assay} · ${esc(color)}</td><td><div aria-hidden="true" class="forest-track forest-${assay.toLowerCase()}">${line}${point}</div></td><td>${Number.isFinite(r)?r.toFixed(3):'—'} ${range}</td><td>${stat(p)}</td><td>${stat(q)}</td><td>${nHD??'—'} / ${nTarget??'—'}</td><td>${esc(status)}</td></tr>`}).join('');
  container.innerHTML=`<p class="sub">Squares show point-biserial r. Lines show approximate 95% Fisher-z intervals. q is BH FDR within each network. Low-support estimates are descriptive.</p><div class="table-wrap overlap-forest"><table aria-label="Module diagnosis forest plot"><thead><tr><th>Contrast</th><th>Module</th><th>Point-biserial r</th><th>r [95% CI]</th><th>p</th><th>q</th><th>n HD / disease</th><th>Support</th></tr></thead><tbody>${body}</tbody><tfoot><tr><td colspan="2"></td><td><div class="forest-axis"><span>${-limit}</span><span>0</span><span>${limit}</span></div></td><td colspan="5"></td></tr></tfoot></table></div>`;
}
function renderCrossOmicsExplorer(){
  renderCrossOmicsHeatmap();
  const entry=crossEntry(crossOverlapState.recordId),pair=crossSelectedPair();if(!entry||!pair)return;updateCrossModuleOptions();document.getElementById('crossOmicsRecord').textContent=`${entry.title||crossOverlapState.recordId} · protein=${entry.proteomics_condition||'NA'} · ${entry.specimen||''}`;
  const rows=crossFilteredRows(),pages=Math.max(1,Math.ceil(rows.length/crossOverlapState.pageSize));crossOverlapState.page=Math.min(Math.max(1,crossOverlapState.page),pages);const start=(crossOverlapState.page-1)*crossOverlapState.pageSize,shown=rows.slice(start,start+crossOverlapState.pageSize),badges=[pair.transcript_best?'RNA best':null,pair.proteomics_best?'protein best':null,pair.reciprocal_best?'reciprocal best':null].filter(Boolean).map(label=>`<span class="cross-best-badge">${label}</span>`).join(' ');
  document.getElementById('crossOmicsSummary').innerHTML=`<span>RNA <b>${esc(pair.transcript_module)}</b></span><span>Protein <b>${esc(pair.proteomics_module)}</b></span><span>weighted Jaccard <b>${crossNumber(pair.weighted_jaccard,4)}</b></span><span>plain Jaccard <b>${crossNumber(pair.jaccard,4)}</b></span><span><b>${crossNumber(pair.shared_genes,0)}</b> shared genes</span><span>RNA genes <b>${crossNumber(pair.transcript_module_genes,0)}</b></span><span>protein genes <b>${crossNumber(pair.proteomics_module_genes,0)}</b></span><span>weighted intersection <b>${crossNumber(pair.weighted_intersection,3)}</b></span><span>weighted union <b>${crossNumber(pair.weighted_union,3)}</b></span>${badges}<span>${rows.length.toLocaleString()} search matches</span>`;renderCrossGeneGroups();renderCrossOmicsDiagnosis();
  document.getElementById('crossOmicsRows').innerHTML=shown.map(row=>`<tr><td>${esc(row[2])}</td><td class="num">${crossNumber(row[3])}</td><td class="num">${crossNumber(row[4])}</td><td class="num">${crossNumber(row[5])}</td><td class="num">${crossNumber(row[6])}</td></tr>`).join('')||'<tr><td colspan="5" class="empty">No genes overlap between this module pair.</td></tr>';
  document.getElementById('crossOmicsPage').textContent=rows.length?`${start+1}-${Math.min(start+shown.length,rows.length)} of ${rows.length.toLocaleString()} · page ${crossOverlapState.page}/${pages}`:'0 overlapping genes';document.getElementById('crossOmicsPrev').disabled=crossOverlapState.page<=1;document.getElementById('crossOmicsNext').disabled=crossOverlapState.page>=pages;
}
async function selectCrossPair(recordId,transcriptModule,proteomicsModule){
  const entry=crossEntry(recordId);if(!entry)return;const pairs=crossPairs(entry);let pair=pairs.find(item=>item.transcript_module===transcriptModule&&item.proteomics_module===proteomicsModule)||pairs[0];if(!pair)return;const token=++crossOverlapState.loadingToken;crossOverlapState.recordId=recordId;crossOverlapState.transcriptModule=pair.transcript_module;crossOverlapState.proteomicsModule=pair.proteomics_module;crossOverlapState.payload=null;crossOverlapState.page=1;updateCrossModuleOptions();document.getElementById('crossOmicsRows').innerHTML='<tr><td colspan="5" class="empty">Loading overlapping genes…</td></tr>';
  const proteinRecordId=`proteomics::${entry.specimen}::${entry.proteomics_condition}`,proteinGenes=window.CROSS_OMICS_INDEX?.proteinModules?.[proteinRecordId]?.[pair.proteomics_module]||null;crossOverlapState.transcriptMembers=null;crossOverlapState.proteomicsMembers=proteinGenes;crossOverlapState.geneGroupError=proteinGenes?null:`Mapped Protein module genes are unavailable for ${proteinRecordId} / ${pair.proteomics_module}.`;renderCrossGeneGroups();
  const transcriptMembersPromise=loadModuleShard(recordId,pair.transcript_module);
  renderCrossOmicsHeatmap();renderCrossOmicsDiagnosis();
  try{const payload=await loadCrossShard(recordId);if(token!==crossOverlapState.loadingToken)return;crossOverlapState.payload=payload;renderCrossOmicsExplorer()}catch(error){if(token!==crossOverlapState.loadingToken)return;document.getElementById('crossOmicsRows').innerHTML=`<tr><td colspan="5" class="empty">${esc(error.message)}</td></tr>`}
  try{const members=await transcriptMembersPromise;if(token!==crossOverlapState.loadingToken)return;crossOverlapState.transcriptMembers=members;renderCrossGeneGroups()}catch(error){if(token!==crossOverlapState.loadingToken)return;crossOverlapState.geneGroupError=error.message;renderCrossGeneGroups()}
}
function openCrossOmicsOverlap(recordId,transcriptModule,proteomicsModule){
  const entry=crossEntry(recordId);if(!entry)return;const pairs=crossPairs(entry),pair=pairs.find(item=>item.transcript_module===transcriptModule&&item.proteomics_module===proteomicsModule)||pairs[0];if(!pair)return;document.getElementById('crossOmicsOverlay').hidden=false;document.body.classList.add('cross-omics-open');crossOverlapState.query='';crossOverlapState.sort='intersection_desc';document.getElementById('crossOmicsSearch').value='';document.getElementById('crossOmicsSort').value='intersection_desc';selectCrossPair(recordId,pair.transcript_module,pair.proteomics_module);
}
function closeCrossOmicsOverlap(){document.getElementById('crossOmicsOverlay').hidden=true;document.body.classList.remove('cross-omics-open')}
window.openCrossOmicsOverlap=openCrossOmicsOverlap;
window.addEventListener('message',event=>{const frame=document.querySelector('.final-jaccard-page'),message=event.data;if(event.source!==frame?.contentWindow||message?.type!=='open-cross-omics-overlap'||!crossEntry(message.recordId))return;openCrossOmicsOverlap(message.recordId,message.transcriptModule||null,message.proteomicsModule||null)});
function crossCsvCell(value){const text=String(value??'');return /[",\n]/.test(text)?`"${text.replaceAll('"','""')}"`:text}
function exportCrossOmicsRows(){const pair=crossSelectedPair(),rows=crossFilteredRows(),head=['transcript_module','proteomics_module','gene','transcript_abs_kME','proteomics_abs_kME','weighted_intersection','weighted_union','kME_weighted_jaccard'];const data=rows.map(row=>[row[0],row[1],row[2],row[3],row[4],row[5],row[6],pair?.weighted_jaccard]);const blob=new Blob(['\ufeff'+[head,...data].map(row=>row.map(crossCsvCell).join(',')).join('\n')],{type:'text/csv;charset=utf-8'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=`${crossOverlapState.recordId.replace(/[^a-zA-Z0-9_-]+/g,'_')}__${crossOverlapState.transcriptModule}__${crossOverlapState.proteomicsModule}__overlap_genes.csv`;a.click();URL.revokeObjectURL(a.href)}
document.getElementById('crossOmicsClose').onclick=closeCrossOmicsOverlap;document.getElementById('crossOmicsOverlay').onclick=event=>{if(event.target.id==='crossOmicsOverlay')closeCrossOmicsOverlap()};document.addEventListener('keydown',event=>{if(event.key==='Escape'&&!document.getElementById('crossOmicsOverlay').hidden)closeCrossOmicsOverlap()});
document.getElementById('crossTranscriptModule').onchange=event=>{const entry=crossEntry(crossOverlapState.recordId),pair=crossPairs(entry).find(item=>item.transcript_module===event.target.value);if(pair)selectCrossPair(crossOverlapState.recordId,pair.transcript_module,pair.proteomics_module)};document.getElementById('crossProteomicsModule').onchange=event=>selectCrossPair(crossOverlapState.recordId,crossOverlapState.transcriptModule,event.target.value);document.getElementById('crossOmicsSearch').oninput=event=>{crossOverlapState.query=event.target.value;crossOverlapState.page=1;renderCrossOmicsExplorer()};document.getElementById('crossOmicsSort').onchange=event=>{crossOverlapState.sort=event.target.value;crossOverlapState.page=1;renderCrossOmicsExplorer()};document.getElementById('crossOmicsPrev').onclick=()=>{crossOverlapState.page--;renderCrossOmicsExplorer()};document.getElementById('crossOmicsNext').onclick=()=>{crossOverlapState.page++;renderCrossOmicsExplorer()};document.getElementById('crossOmicsExport').onclick=exportCrossOmicsRows;
document.getElementById('crossOmicsHeatmap').onclick=event=>{const cell=event.target.closest('.cross-heatmap-cell');if(cell)selectCrossPair(crossOverlapState.recordId,cell.dataset.crossRna,cell.dataset.crossProtein)};
document.addEventListener('click',event=>{const button=event.target.closest('[data-open-cross-overlap]');if(button)openCrossOmicsOverlap(button.dataset.openCrossOverlap,button.dataset.transcriptModule||null,button.dataset.proteomicsModule||null)});
let crossDecoratePending=false;
function decorateCrossOmicsEntrypoints(){
  crossDecoratePending=false;
  if(['cross_omics_similarity','after_outlier_cross_omics_similarity'].includes(comparisonArtifactChoice))document.querySelectorAll('.condition-panel').forEach(panel=>{const inspect=panel.querySelector('[data-open-record]'),id=inspect?.dataset.openRecord,variantId=crossVariantRecordId(id,comparisonArtifactChoice);if(!id||!crossEntry(variantId)||panel.querySelector('[data-open-cross-overlap]'))return;const button=document.createElement('button');button.type='button';button.className='cross-overlap-button';button.dataset.openCrossOverlap=variantId;button.textContent='View shared genes';panel.querySelector('.condition-links')?.appendChild(button)});
  const detailActions=document.querySelector('#detail .detail-head .actions'),detailArtifact=artifactChoice[selectedId],detailVariant=crossVariantRecordId(selectedId,detailArtifact);if(detailActions&&['cross_omics_similarity','after_outlier_cross_omics_similarity'].includes(detailArtifact)&&crossEntry(detailVariant)&&!detailActions.querySelector('[data-open-cross-overlap]')){const button=document.createElement('button');button.type='button';button.dataset.openCrossOverlap=detailVariant;button.textContent='View shared genes';detailActions.appendChild(button)}
  document.querySelectorAll('#compareGrid .compare-card').forEach((card,index)=>{const id=compareIds[index];if(!id||!crossEntry(id)||card.querySelector('[data-open-cross-overlap]'))return;const button=document.createElement('button');button.type='button';button.className='cross-overlap-button';button.dataset.openCrossOverlap=id;button.textContent='View shared genes';card.appendChild(button)});
}
function scheduleCrossDecoration(){if(crossDecoratePending)return;crossDecoratePending=true;queueMicrotask(decorateCrossOmicsEntrypoints)}
new MutationObserver(scheduleCrossDecoration).observe(document.getElementById('cellComparison').parentElement,{childList:true,subtree:true});scheduleCrossDecoration();
const crossParams=new URLSearchParams(location.search),deepCrossRecord=crossParams.get('cross_record');if(deepCrossRecord&&crossEntry(deepCrossRecord))setTimeout(()=>openCrossOmicsOverlap(deepCrossRecord,crossParams.get('transcript_module'),crossParams.get('proteomics_module')),0);
'''


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
  <select id="viewMode"><option value="cell">Cell type comparison</option><option value="outliers">Outlier candidates</option><option value="run">Individual run inspection</option></select><select id="source"><option value="">All sources</option></select><select id="cell"><option value="">Choose cell type</option></select><select class="run-only" id="threshold"><option value="">All thresholds</option></select><select class="run-only" id="strategy"><option value="">All strategies</option></select><select class="run-only" id="mode"><option value="">All modes</option></select><select class="run-only" id="specimen"><option value="">All specimens</option></select><select class="run-only" id="status"><option value="">All statuses</option></select><select class="run-only" id="decision"><option value="">All decisions</option><option value="unreviewed">Unreviewed</option><option value="use">Use</option><option value="hold">Hold</option><option value="reject">Do not use</option></select><select class="run-only" id="sort"><option value="default">Default order</option><option value="modules_desc">Final modules ↓</option><option value="unassigned_desc">Unassigned ↓</option><option value="runtime_desc">Runtime ↓</option></select><input class="run-only" id="search" placeholder="Search cell type or network"><span class="spacer"></span><button class="run-only" id="nextUnreviewed">Next unreviewed</button><button class="run-only" id="showCompare">Compare <span id="compareCount">0</span></button><button id="exportReview">Export review CSV</button><label><button type="button" id="importButton">Import review CSV</button><input id="importReview" type="file" accept=".csv" hidden></label>
</div>
<main class="cell-comparison cell-only" id="cellComparison"></main>
<main class="outlier-page outlier-only" id="outlierReview"></main>
<section class="compare run-only" id="comparePanel"><div class="detail-head"><h2>Comparison workspace</h2><div class="actions"><button id="clearCompare">Clear</button></div></div><div class="compare-grid" id="compareGrid"></div></section>
<main class="layout run-only" id="runInspection"><aside class="list"><div class="list-head" id="listCount"></div><div id="recordList"></div></aside><section class="detail" id="detail"></section></main>
<script src="dashboard-data.js"></script>
<script>
const bundle=window.EXPERIMENT_DASHBOARD_DATA;const records=bundle.records;const storageKey=`experiment-review:${bundle.schema_version}:${bundle.run_root}`;let reviews=JSON.parse(localStorage.getItem(storageKey)||'{}');let selectedId=records[0]?.id||null;let artifactChoice={};let comparisonArtifactChoice=null;let comparisonHierarchyChoice={specimen:null,threshold:null};let compareIds=[];
const $=id=>document.getElementById(id);const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));const fmt=(v,d=2)=>v===null||v===undefined||v===''?'—':typeof v==='number'?v.toLocaleString(undefined,{maximumFractionDigits:d}):esc(v);const stat=v=>v===null||v===undefined||v===''?'—':typeof v==='number'&&Math.abs(v)>0&&Math.abs(v)<0.0001?v.toExponential(2):fmt(v,4);const pct=v=>v===null||v===undefined?'—':`${(v*100).toFixed(1)}%`;
const finalJaccardOption=document.createElement('option');finalJaccardOption.value='final_jaccard';finalJaccardOption.textContent='Final CSF Jaccard';$('viewMode').insertBefore(finalJaccardOption,$('viewMode').options[1]||null);
const finalJaccardMain=document.createElement('main');finalJaccardMain.className='final-jaccard-only';finalJaccardMain.id='finalJaccardPage';finalJaccardMain.innerHTML='<iframe class="final-jaccard-page" src="final-jaccard.html" title="Final CSF cross-omics Jaccard comparison"></iframe>';document.getElementById('cellComparison').insertAdjacentElement('afterend',finalJaccardMain);
const artifactOrder=['feature_selection','before_heatmap','before_pca','before_pca_3d','after_outlier_heatmap','after_outlier_pca','after_outlier_pca_3d','cut_diagnostic','selected_cut','after_outlier_cut_diagnostic','after_heatmap','after_pca','after_pca_3d','soft_threshold','gene_dendrogram','me_dendrogram','module_sizes','after_outlier_soft_threshold','after_outlier_gene_dendrogram','after_outlier_me_dendrogram','after_outlier_module_sizes','diagnosis_pairwise_hd','diagnosis_heatmap','cross_omics_similarity','after_outlier_diagnosis_pairwise_hd','after_outlier_diagnosis_heatmap','after_outlier_cross_omics_similarity'];
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
function prettyStrategy(v){return ({'01_cpm_detection':'CPM detection','02_hvg_top9000':'HVG top 9,000','03_hvg_top5000':'HVG top 5,000'})[v]||prettyDimension(v)}
function comparisonRows(){
  const source=$('source').value,cell=$('cell').value;
  return records.filter(r=>(!source||r.source_id===source)&&cell&&r.dimensions.cell_type===cell).sort((a,b)=>{
    const numberIn=v=>Number(String(v??'').match(/\d+/)?.[0]??999),specimenRank=v=>v==='CSF'?0:v==='PB'?1:2;
    const ka=[specimenRank(a.dimensions.specimen),numberIn(a.dimensions.threshold),numberIn(a.dimensions.strategy),String(a.dimensions.mode??''),String(a.provenance.condition??'')];
    const kb=[specimenRank(b.dimensions.specimen),numberIn(b.dimensions.threshold),numberIn(b.dimensions.strategy),String(b.dimensions.mode??''),String(b.provenance.condition??'')];
    for(let i=0;i<ka.length;i++){const d=typeof ka[i]==='number'?ka[i]-kb[i]:ka[i].localeCompare(kb[i]);if(d)return d}return 0;
  });
}
function comparisonHierarchy(rows){
  const ordered=(values,key)=>[...new Set(values.filter(Boolean))].sort((a,b)=>{
    if(key==='specimen'){const rank=v=>v==='CSF'?0:v==='PB'?1:2;return rank(a)-rank(b)||String(a).localeCompare(String(b))}
    const number=v=>Number(String(v).match(/\d+/)?.[0]??999);return number(a)-number(b)||String(a).localeCompare(String(b));
  });
  const specimens=ordered(rows.map(r=>r.dimensions.specimen),'specimen');
  if(!specimens.includes(comparisonHierarchyChoice.specimen))comparisonHierarchyChoice.specimen=specimens[0]||null;
  const specimenRows=rows.filter(r=>!comparisonHierarchyChoice.specimen||r.dimensions.specimen===comparisonHierarchyChoice.specimen);
  const thresholds=ordered(specimenRows.map(r=>r.dimensions.threshold),'threshold');
  if(!thresholds.includes(comparisonHierarchyChoice.threshold))comparisonHierarchyChoice.threshold=thresholds[0]||null;
  const thresholdRows=specimenRows.filter(r=>!comparisonHierarchyChoice.threshold||r.dimensions.threshold===comparisonHierarchyChoice.threshold);
  return {specimens,thresholds,shown:thresholdRows};
}
function comparisonDimension(label,key,values){
  if(!values.length)return '';
  return `<div class="comparison-dimension"><div class="comparison-dimension-label">${esc(label)}</div><div class="comparison-dimension-options">${values.map(value=>`<button type="button" data-comparison-${esc(key)}="${esc(value)}" class="${comparisonHierarchyChoice[key]===value?'active':''}">${esc(prettyDimension(value))}</button>`).join('')}</div></div>`;
}
const crossOmicsArtifactIds=new Set(['cross_omics_similarity','after_outlier_cross_omics_similarity']);
function isOtherCellRecord(record){return String(record?.dimensions?.cell_type||'').trim().toLowerCase()==='other'}
function visibleArtifacts(record){return record.artifacts.filter(artifact=>!(isOtherCellRecord(record)&&crossOmicsArtifactIds.has(artifact.id)))}
for(const record of records)if(isOtherCellRecord(record))record.artifacts=visibleArtifacts(record);
function comparisonArtifacts(rows){
  const byId=new Map();for(const r of rows)for(const a of visibleArtifacts(r))if(!byId.has(a.id))byId.set(a.id,a);
  return [...byId.values()].sort((a,b)=>{const ai=artifactOrder.indexOf(a.id),bi=artifactOrder.indexOf(b.id);return (ai<0?999:ai)-(bi<0?999:bi)||a.label.localeCompare(b.label)});
}
function renderCellComparison(){
  const rows=comparisonRows(),artifacts=comparisonArtifacts(rows),cell=$('cell').value;
  if(!cell){$('cellComparison').innerHTML='<div class="empty">Choose one cell type to compare its conditions.</div>';return}
  if(!rows.length){$('cellComparison').innerHTML='<div class="empty">No experiment conditions are available for this source and cell type.</div>';return}
  const hierarchy=comparisonHierarchy(rows),shownRows=hierarchy.shown;
  if(!comparisonArtifactChoice||!artifacts.some(a=>a.id===comparisonArtifactChoice))comparisonArtifactChoice=artifacts[0]?.id||null;
  const chosenMeta=artifacts.find(a=>a.id===comparisonArtifactChoice),available=shownRows.filter(r=>visibleArtifacts(r).some(a=>a.id===comparisonArtifactChoice)).length;
  const tabs=artifacts.map(a=>`<button data-comparison-artifact="${esc(a.id)}" class="${a.id===comparisonArtifactChoice?'active':''}">${esc(a.label)}</button>`).join('');
  const dimensionNav=[comparisonDimension('Specimen','specimen',hierarchy.specimens),comparisonDimension('Cell threshold','threshold',hierarchy.thresholds)].join('');
  const panels=shownRows.map(r=>{
    const a=visibleArtifacts(r).find(x=>x.id===comparisonArtifactChoice),rv=reviewOf(r.id);
    const afterOutlier=String(comparisonArtifactChoice||'').startsWith('after_outlier_');
    const shownSamples=afterOutlier?(r.metrics.after_outlier_samples??r.metrics.samples):r.metrics.samples;
    const shownGenes=afterOutlier?(r.metrics.after_outlier_genes??r.metrics.genes):r.metrics.genes;
    const shownDynamic=afterOutlier?(r.metrics.after_outlier_dynamic_modules??r.metrics.dynamic_modules):r.metrics.dynamic_modules;
    const shownFinal=afterOutlier?(r.metrics.after_outlier_final_modules??r.metrics.final_modules):r.metrics.final_modules;
    const condition=prettyStrategy(r.dimensions.strategy);
    const secondary=[r.dimensions.mode,`n=${fmt(shownSamples)}`,`genes=${fmt(shownGenes)}`,`modules=${fmt(shownDynamic)}→${fmt(shownFinal)}`].filter(Boolean).map(prettyDimension).join(' · ');
    const visual=!a?'<div class="condition-missing">This artifact was not generated for this condition.</div>':a.kind==='html'?`<iframe loading="lazy" src="${esc(a.path)}" title="${esc(a.label)} · ${esc(condition)}"></iframe>`:`<a href="${esc(a.path)}" target="_blank"><img loading="lazy" src="${esc(a.path)}" alt="${esc(a.label)} · ${esc(condition)}"></a>`;
    const open=a?`<a href="${esc(a.path)}" target="_blank">Open artifact</a>`:'';
    return `<article class="condition-panel"><div class="condition-head"><div><h3>${esc(condition||r.title)}</h3><div class="condition-meta">${secondary}</div></div><span class="condition-review ${esc(rv.decision)}">${esc(rv.decision)}</span></div><div class="condition-viewer">${visual}</div><div class="condition-links">${open}<button data-open-record="${esc(r.id)}">Inspect this run</button></div></article>`;
  }).join('');
  const path=[comparisonHierarchyChoice.specimen,comparisonHierarchyChoice.threshold].filter(Boolean).map(prettyDimension).join(' → ');
  $('cellComparison').innerHTML=`<div class="comparison-head"><div><h2>${esc(prettyDimension(cell))}</h2><div class="sub">Choose specimen and cell threshold, then compare all feature strategies side by side.</div></div><div class="comparison-counts">${shownRows.length} of ${rows.length} conditions · ${artifacts.length} output tabs · ${available}/${shownRows.length} shown</div></div><div class="comparison-hierarchy" aria-label="Experiment condition hierarchy">${dimensionNav}</div><div class="artifact-tabs">${tabs}</div><div class="comparison-note">Current condition: <b>${esc(path||'—')}</b> · feature strategies are aligned in parallel · current output: <b>${esc(chosenMeta?.label||'—')}</b></div><section class="condition-grid strategy-grid">${panels}</section>`;
  document.querySelectorAll('[data-comparison-specimen]').forEach(n=>n.onclick=()=>{comparisonHierarchyChoice.specimen=n.dataset.comparisonSpecimen;comparisonHierarchyChoice.threshold=null;renderCellComparison()});
  document.querySelectorAll('[data-comparison-threshold]').forEach(n=>n.onclick=()=>{comparisonHierarchyChoice.threshold=n.dataset.comparisonThreshold;renderCellComparison()});
  document.querySelectorAll('[data-comparison-artifact]').forEach(n=>n.onclick=()=>{comparisonArtifactChoice=n.dataset.comparisonArtifact;renderCellComparison()});
  document.querySelectorAll('[data-open-record]').forEach(n=>n.onclick=()=>{selectedId=n.dataset.openRecord;$('viewMode').value='run';renderAll()});
}
function renderKpis(){const complete=records.filter(r=>r.status==='complete_findModules').length,collapsed=records.filter(r=>r.metrics.collapsed).length,reviewed=records.filter(r=>reviewOf(r.id).decision!=='unreviewed').length,holds=records.filter(r=>reviewOf(r.id).decision==='hold').length,uses=records.filter(r=>reviewOf(r.id).decision==='use').length;const values=[[records.length,'Experiment records'],[complete,'findModules complete'],[collapsed,'Collapsed'],[reviewed,'Reviewed'],[uses,'Use'],[holds,'Hold']];$('kpis').innerHTML=values.map(([v,l])=>`<div class="kpi"><b>${v}</b><span>${l}</span></div>`).join('')}
function renderList(){const rows=filtered();$('listCount').textContent=`${rows.length} / ${records.length} records`;if(rows.length&&!rows.some(r=>r.id===selectedId))selectedId=rows[0].id;$('recordList').innerHTML=rows.map(r=>{const rv=reviewOf(r.id),warn=r.flags.some(f=>f.level==='warning');return `<button class="record ${r.id===selectedId?'active':''}" data-id="${esc(r.id)}"><div class="record-title">${esc(r.title)} ${warn?'<span class="pill warn">flag</span>':''}</div><div class="record-meta"><span>${esc(r.dimensions.threshold)}</span><span>${esc(r.dimensions.strategy)}</span><span>${esc(r.dimensions.mode)}</span></div><div class="record-metrics"><span>n=${fmt(r.metrics.samples)}</span><span>genes=${fmt(r.metrics.genes)}</span><span>dynamic=${fmt(r.metrics.dynamic_modules)}</span><span>final=${fmt(r.metrics.final_modules)}</span><span class="pill ${rv.decision==='use'?'ok':rv.decision==='hold'||rv.decision==='reject'?'warn':''}">${esc(rv.decision)}</span></div></button>`}).join('')||'<div class="empty">No matching experiments</div>';document.querySelectorAll('.record').forEach(n=>n.onclick=()=>{selectedId=n.dataset.id;renderAll()})}
function metric(label,value){return `<div class="metric"><div class="label">${esc(label)}</div><div class="value">${value}</div></div>`}
function renderDetail(){const r=records.find(x=>x.id===selectedId);if(!r){$('detail').innerHTML='<div class="empty">Select an experiment</div>';return}const rv=reviewOf(r.id),artifacts=r.artifacts;if(!artifactChoice[r.id]||!artifacts.some(a=>a.id===artifactChoice[r.id]))artifactChoice[r.id]=(artifacts.find(a=>a.id==='diagnosis_heatmap')||artifacts.find(a=>a.id==='me_dendrogram')||artifacts[0])?.id;const chosen=artifacts.find(a=>a.id===artifactChoice[r.id]);const metrics=[['Samples',fmt(r.metrics.samples)],['Genes',fmt(r.metrics.genes)],['Power',fmt(r.metrics.selected_power)],['Selected SFT R²',fmt(r.metrics.selected_power_sft_r2,3)],['Mean connectivity',fmt(r.metrics.selected_power_mean_connectivity,1)],['Dynamic modules',fmt(r.metrics.dynamic_modules)],['Final modules',fmt(r.metrics.final_modules)],['Unassigned',pct(r.metrics.unassigned_fraction)],['Largest module',pct(r.metrics.largest_module_fraction)],['Median |kME|',fmt(r.metrics.median_abs_kme,3)],['Diagnosis tests',fmt(r.metrics.diagnosis_tests)],['Min diagnosis q',stat(r.metrics.minimum_diagnosis_q)],['Samples cut',fmt(r.metrics.samples_removed_by_cut)],['Runtime (s)',fmt(r.metrics.runtime_seconds,1)]];const groups={};for(const f of r.files)(groups[f.category]??=[]).push(f);const proper=r.modules.filter(m=>!m.fully_unassigned).sort((a,b)=>(b.features??0)-(a.features??0));$('detail').innerHTML=`<div class="detail-head"><div><h2>${esc(r.title)}</h2><div class="sub">${esc(r.provenance.network_id)} · ${esc(r.status)}</div></div><div class="actions"><button id="addCompare">${compareIds.includes(r.id)?'Remove from compare':'Add to compare'}</button></div></div><div class="flags">${r.flags.map(f=>`<span class="flag">${esc(f.label)}</span>`).join('')}</div><div class="metric-grid">${metrics.map(x=>metric(...x)).join('')}</div><section class="review"><div><b>Decision</b><div class="decision-row">${[['use','Use'],['hold','Hold'],['reject','Do not use'],['unreviewed','Clear']].map(([v,l])=>`<button data-decision="${v}" class="${rv.decision===v?'selected':''}">${l}</button>`).join('')}</div><div class="checklist"><label><input type="checkbox" data-check="power" ${rv.checks.power?'checked':''}> Soft-threshold behavior inspected</label><label><input type="checkbox" data-check="modules" ${rv.checks.modules?'checked':''}> ME merging and module sizes inspected</label><label><input type="checkbox" data-check="outliers" ${rv.checks.outliers?'checked':''}> Sample/outlier evidence inspected</label><label><input type="checkbox" data-check="ready" ${rv.checks.ready?'checked':''}> Ready for downstream analysis</label></div></div><div><b>Review notes</b><textarea id="reviewNotes" placeholder="Why use, hold, or reject this run?">${esc(rv.notes)}</textarea><div class="sub">Saved in this browser · ${rv.updated_at?`updated ${esc(rv.updated_at)}`:'not reviewed yet'}</div></div></section><h3 class="section-title">Visual evidence</h3><div class="artifact-tabs">${artifacts.map(a=>`<button data-artifact="${a.id}" class="${a.id===chosen?.id?'active':''}">${esc(a.label)}</button>`).join('')}</div><div class="viewer">${!chosen?'<div class="empty">No visual artifact</div>':chosen.kind==='html'?`<div style="width:100%"><div class="open"><a href="${esc(chosen.path)}" target="_blank">Open interactive result in a new tab</a></div><iframe src="${esc(chosen.path)}"></iframe></div>`:`<a href="${esc(chosen.path)}" target="_blank"><img src="${esc(chosen.path)}" alt="${esc(chosen.label)}"></a>`}</div><div class="lower"><div><h3 class="section-title">Final proper modules</h3><div class="table-wrap"><table><thead><tr><th>Module</th><th>Features</th><th>Median |kME|</th></tr></thead><tbody>${proper.map(m=>`<tr><td><span class="swatch" style="background:${esc(m.color)}"></span>${esc(m.color)}</td><td>${fmt(m.features)}</td><td>${fmt(m.median_abs_kme,3)}</td></tr>`).join('')}</tbody></table></div></div><div><h3 class="section-title">Files and provenance</h3><div class="files">${Object.entries(groups).map(([g,fs])=>`<details ${g==='summary'?'open':''}><summary>${esc(g)} (${fs.length})</summary>${fs.map(f=>`<a href="${esc(f.path)}" target="_blank">${esc(f.label)}</a>`).join('')}</details>`).join('')}<details><summary>Parameters</summary>${Object.entries(r.provenance).map(([k,v])=>`<div><b>${esc(k)}</b>: ${fmt(v)}</div>`).join('')}</details></div></div></div>`;document.querySelectorAll('[data-decision]').forEach(n=>n.onclick=()=>saveReview(r.id,{decision:n.dataset.decision}));document.querySelectorAll('[data-check]').forEach(n=>n.onchange=()=>saveReview(r.id,{checks:{...reviewOf(r.id).checks,[n.dataset.check]:n.checked}}));$('reviewNotes').onchange=e=>saveReview(r.id,{notes:e.target.value});document.querySelectorAll('[data-artifact]').forEach(n=>n.onclick=()=>{artifactChoice[r.id]=n.dataset.artifact;renderDetail()});$('addCompare').onclick=()=>{compareIds=compareIds.includes(r.id)?compareIds.filter(x=>x!==r.id):[...compareIds,r.id].slice(-4);renderAll()}}
function renderCompare(){const panel=$('comparePanel');$('compareCount').textContent=compareIds.length;if(!compareIds.length){panel.classList.remove('visible');$('compareGrid').innerHTML='';return}const selected=compareIds.map(id=>records.find(r=>r.id===id)).filter(Boolean);$('compareGrid').innerHTML=selected.map(r=>{const a=r.artifacts.find(x=>x.id==='me_dendrogram')||r.artifacts[0];return `<div class="compare-card"><h3>${esc(r.title)}</h3><div class="sub">${esc(r.dimensions.strategy)} · dynamic ${fmt(r.metrics.dynamic_modules)} → final ${fmt(r.metrics.final_modules)}</div>${a?`<a href="${esc(a.path)}" target="_blank"><img src="${esc(a.path)}" alt="${esc(a.label)}"></a>`:''}</div>`}).join('')}
function renderAll(){const view=$('viewMode').value;document.body.dataset.view=view;renderKpis();if(view==='cell'){renderCellComparison()}else if(view==='outliers'){window.OutlierReview?.render($('outlierReview'));window.OutlierReview?.renderKpis()}else if(view!=='final_jaccard'){renderList();renderDetail();renderCompare();window.OutlierReview?.attachRunShortcut(selectedId)}$('sourceLine').textContent=`Schema ${bundle.schema_version} · generated ${bundle.generated_at_utc} · ${bundle.sources.map(s=>s.label).join(', ')}`}
window.renderAll=renderAll;
window.openOutlierRecord=(recordId,artifactId)=>{selectedId=recordId;if(artifactId)artifactChoice[recordId]=artifactId;$('viewMode').value='run';renderAll()};
$('showCompare').onclick=()=>{$('comparePanel').classList.toggle('visible');$('comparePanel').scrollIntoView({behavior:'smooth'})};$('clearCompare').onclick=()=>{compareIds=[];renderAll()};$('nextUnreviewed').onclick=()=>{const rows=filtered(),i=Math.max(-1,rows.findIndex(r=>r.id===selectedId));const next=[...rows.slice(i+1),...rows.slice(0,i+1)].find(r=>reviewOf(r.id).decision==='unreviewed');if(next){selectedId=next.id;renderAll()}};
function csvCell(v){const s=String(v??'');return /[",\n]/.test(s)?`"${s.replaceAll('"','""')}"`:s}function exportReview(){const head=['record_id','source_id','network_id','decision','power_checked','modules_checked','outliers_checked','ready_checked','notes','updated_at'];const rows=records.map(r=>{const v=reviewOf(r.id);return [r.id,r.source_id,r.provenance.network_id,v.decision,v.checks.power,v.checks.modules,v.checks.outliers,v.checks.ready,v.notes,v.updated_at]}).filter(x=>x[3]!=='unreviewed'||x[8]).map(x=>x.map(csvCell).join(','));const blob=new Blob([[head.join(','),...rows].join('\n')],{type:'text/csv;charset=utf-8'});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=`experiment_review_${new Date().toISOString().slice(0,10)}.csv`;a.click();URL.revokeObjectURL(a.href)}$('exportReview').onclick=exportReview;$('importButton').onclick=()=>$('importReview').click();
function parseCSV(text){const rows=[];let row=[],cell='',quoted=false;for(let i=0;i<text.length;i++){const c=text[i];if(quoted&&c==='"'&&text[i+1]==='"'){cell+='"';i++}else if(c==='"'){quoted=!quoted}else if(c===','&&!quoted){row.push(cell);cell=''}else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(cell);if(row.some(x=>x!==''))rows.push(row);row=[];cell=''}else cell+=c}row.push(cell);if(row.some(x=>x!==''))rows.push(row);return rows}$('importReview').onchange=e=>{const file=e.target.files[0];if(!file)return;const reader=new FileReader();reader.onload=()=>{const rows=parseCSV(reader.result),head=rows.shift(),ix=Object.fromEntries(head.map((x,i)=>[x,i]));for(const r of rows){const id=r[ix.record_id];if(!records.some(x=>x.id===id))continue;reviews[id]={decision:r[ix.decision]||'unreviewed',checks:{power:r[ix.power_checked]==='true',modules:r[ix.modules_checked]==='true',outliers:r[ix.outliers_checked]==='true',ready:r[ix.ready_checked]==='true'},notes:r[ix.notes]||'',updated_at:r[ix.updated_at]||new Date().toISOString()}}localStorage.setItem(storageKey,JSON.stringify(reviews));renderAll()};reader.readAsText(file,'utf-8')};const requestedView=new URLSearchParams(location.search).get('view');if([...$('viewMode').options].some(option=>option.value===requestedView))$('viewMode').value=requestedView;renderAll();
</script>
</body>
</html>'''


if __name__ == "__main__":
    arguments = parse_args()
    build(arguments.run_root.resolve())
