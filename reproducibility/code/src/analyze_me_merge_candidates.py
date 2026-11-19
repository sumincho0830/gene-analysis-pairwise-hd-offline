"""Diagnose transcriptomics module merges before rerunning findModules.

This script does not rerun PyWGCNA.  It reuses the saved pre-merge dynamic
module eigengenes and their average-linkage tree to compare the baseline
MEDissThres=0.20 with a stricter diagnostic cut at 0.10.  A current final
module is a candidate when its constituent dynamic modules occupy more than
one component at the 0.10 cut.

The diagnostic components are not claimed to be the exact final modules from
a new PyWGCNA run: mergeCloseModules is iterative and recomputes eigengenes
after each merge.  Exact module assignments require a separate rerun.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster

from pipeline_provenance import write_run_manifest


BASELINE_THRESHOLD = 0.20
CANDIDATE_THRESHOLD = 0.10
CONDITION = "01_default_ds2_mms50_med0p20"
ACTUAL_CONDITION = "02_predicted_ds2_mms50_med0p10"
OUTPUT_DIRECTORY = "07_module_merge_diagnostics"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def clean_text(value: Any, fallback: str = "ALL") -> str:
    if pd.isna(value) or str(value).strip() == "":
        return fallback
    return str(value)


def network_root(stage_root: Path, row: pd.Series) -> Path:
    root = (
        stage_root
        / clean_text(row["threshold"])
        / clean_text(row["strategy"])
        / clean_text(row["mode"])
        / clean_text(row["specimen"])
    )
    if clean_text(row["mode"]) == "cell_separate":
        root = root / clean_text(row["cell_type"])
    return root


def network_key(row: pd.Series) -> str:
    return "|".join(
        clean_text(row[column])
        for column in ("threshold", "strategy", "mode", "specimen", "cell_type")
    )


def safe_filename(value: str) -> str:
    readable = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")[:90]
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:10]
    return f"{readable}_{digest}.png"


def relative_link(path: Path, origin: Path) -> str:
    return os.path.relpath(path, start=origin).replace("\\", "/")


def priority_label(
    *,
    current_final_count: int,
    potential_additional: int,
    component_count: int,
) -> str:
    if current_final_count <= 1 and potential_additional > 0:
        return "high"
    if potential_additional >= 3 or component_count >= 4:
        return "high"
    if potential_additional >= 2:
        return "medium"
    return "low"


def pairwise_dissimilarity_stats(
    dynamic_mes: pd.DataFrame, dynamic_columns: list[str]
) -> tuple[float | None, float | None, float | None]:
    if len(dynamic_columns) < 2:
        return None, None, None
    correlations = dynamic_mes.loc[:, dynamic_columns].corr().to_numpy(dtype=float)
    dissimilarities = 1.0 - correlations
    values = dissimilarities[np.triu_indices_from(dissimilarities, k=1)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None, None, None
    return float(values.min()), float(np.median(values)), float(values.max())


def diagnostic_components(
    colors: list[str], cluster_by_color: dict[str, int], sizes: dict[str, int]
) -> list[dict[str, Any]]:
    groups: dict[int, list[str]] = {}
    for color in colors:
        groups.setdefault(int(cluster_by_color[color]), []).append(color)
    return [
        {
            "diagnostic_cluster": int(cluster),
            "dynamic_modules": members,
            "feature_count": int(sum(sizes[color] for color in members)),
        }
        for cluster, members in sorted(groups.items())
    ]


def draw_diagnostic_figure(
    *,
    linkage_matrix: np.ndarray,
    labels: list[str],
    figure_path: Path,
    title: str,
    current_final_count: int,
    candidate_module_count: int,
    potential_additional: int,
) -> None:
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    width = max(11.0, min(22.0, 8.0 + len(labels) * 0.34))
    figure, axis = plt.subplots(figsize=(width, 7.2), facecolor="white")
    if len(labels) == 1:
        axis.set_xlim(-0.75, 0.75)
        axis.set_xticks([0])
        axis.set_xticklabels([f"ME{labels[0]}"])
        axis.text(
            0,
            0.145,
            "Only one proper dynamic module; no ME merge can be split",
            ha="center",
            va="center",
            color="#666666",
        )
    else:
        dendrogram(
            linkage_matrix,
            labels=[f"ME{label}" for label in labels],
            color_threshold=CANDIDATE_THRESHOLD,
            above_threshold_color="#8a8a8a",
            leaf_rotation=90,
            leaf_font_size=8,
            ax=axis,
        )
    axis.axhline(
        BASELINE_THRESHOLD,
        color="#262626",
        linewidth=1.25,
        linestyle="--",
        label="current MEDissThres = 0.20",
    )
    axis.axhline(
        CANDIDATE_THRESHOLD,
        color="#007c83",
        linewidth=1.5,
        linestyle="-",
        label="candidate MEDissThres = 0.10",
    )
    for label in axis.get_xticklabels():
        color = label.get_text()[2:]
        if mcolors.is_color_like(color):
            label.set_color(color)
    axis.set_title(title, loc="left", fontsize=13, fontweight="normal")
    axis.set_ylabel("ME dissimilarity (1 - Pearson correlation)")
    axis.set_xlabel("Dynamic modules before ME merging")
    lower, upper = axis.get_ylim()
    axis.set_ylim(lower, max(upper, BASELINE_THRESHOLD * 1.12))
    axis.text(
        0.995,
        0.98,
        (
            f"current final modules: {current_final_count}\n"
            f"candidate final modules: {candidate_module_count}\n"
            f"potential additional modules: {potential_additional}"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#555555",
        bbox={"facecolor": "white", "edgecolor": "#dddddd", "pad": 5},
    )
    axis.legend(loc="upper left", frameon=False, fontsize=9)
    figure.tight_layout()
    figure.savefig(figure_path, dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def analyze_network(
    *,
    stage_root: Path,
    output_root: Path,
    row: pd.Series,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[Path]]:
    key = network_key(row)
    root = network_root(stage_root, row)
    module_root = root / "05_find_modules" / CONDITION
    analysis_path = module_root / "analysis_summary.json"
    assignments_path = module_root / "module_assignments.csv.gz"
    dynamic_me_path = module_root / "dynamic_module_eigengenes.csv"
    linkage_path = module_root / "module_eigengene_linkage.csv"
    required = (analysis_path, assignments_path, dynamic_me_path, linkage_path)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"{key}: missing required files: {missing}")

    analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
    assignments = pd.read_csv(assignments_path)
    dynamic_mes = pd.read_csv(dynamic_me_path, index_col="sample_id")
    linkage_table = pd.read_csv(linkage_path)

    linkage_matrix = linkage_table[
        ["left_cluster", "right_cluster", "distance", "leaf_count"]
    ].to_numpy(dtype=float)
    if linkage_matrix.shape != (dynamic_mes.shape[1] - 1, 4):
        raise ValueError(f"{key}: linkage/ME dimension mismatch")

    labels = [column[2:] if column.startswith("ME") else column for column in dynamic_mes]
    if len(set(labels)) != len(labels):
        raise ValueError(f"{key}: duplicate dynamic module labels")
    proper = assignments.loc[~assignments["is_raw_unassigned"].astype(bool)].copy()
    assignment_colors = set(proper["dynamic_color"].astype(str))
    if assignment_colors != set(labels):
        raise ValueError(f"{key}: dynamic ME labels do not match assignments")

    if len(labels) == 1:
        clusters_020 = np.ones(1, dtype=int)
        clusters_010 = np.ones(1, dtype=int)
    else:
        clusters_020 = fcluster(
            linkage_matrix, t=BASELINE_THRESHOLD, criterion="distance"
        ).astype(int)
        clusters_010 = fcluster(
            linkage_matrix, t=CANDIDATE_THRESHOLD, criterion="distance"
        ).astype(int)
    cluster020_by_color = dict(zip(labels, clusters_020, strict=True))
    cluster010_by_color = dict(zip(labels, clusters_010, strict=True))
    dynamic_sizes = (
        proper.groupby(proper["dynamic_color"].astype(str), sort=False)
        .size()
        .astype(int)
        .to_dict()
    )

    current_final_count = int(analysis["final_module_count"])
    module_rows: list[dict[str, Any]] = []
    for final_color, group in proper.groupby("module_color", sort=False):
        dynamic_colors = list(dict.fromkeys(group["dynamic_color"].astype(str)))
        components_010 = diagnostic_components(
            dynamic_colors, cluster010_by_color, dynamic_sizes
        )
        components_020 = diagnostic_components(
            dynamic_colors, cluster020_by_color, dynamic_sizes
        )
        component_count = len(components_010)
        potential_additional = component_count - 1
        minimum, median, maximum = pairwise_dissimilarity_stats(
            dynamic_mes, [f"ME{color}" for color in dynamic_colors]
        )
        module_rows.append(
            {
                "network_key": key,
                "threshold": clean_text(row["threshold"]),
                "strategy": clean_text(row["strategy"]),
                "mode": clean_text(row["mode"]),
                "specimen": clean_text(row["specimen"]),
                "cell_type": clean_text(row["cell_type"]),
                "current_final_module": str(final_color),
                "feature_count": int(len(group)),
                "dynamic_module_count": int(len(dynamic_colors)),
                "dynamic_modules": ";".join(dynamic_colors),
                "dynamic_module_sizes": ";".join(
                    f"{color}:{dynamic_sizes[color]}" for color in dynamic_colors
                ),
                "diagnostic_components_at_0p20": int(len(components_020)),
                "diagnostic_components_at_0p10": int(component_count),
                "potential_additional_modules": int(potential_additional),
                "pairwise_me_dissimilarity_min": minimum,
                "pairwise_me_dissimilarity_median": median,
                "pairwise_me_dissimilarity_max": maximum,
                "component_detail_json": json.dumps(
                    components_010, ensure_ascii=False, separators=(",", ":")
                ),
                "candidate": bool(potential_additional > 0),
                "priority": (
                    priority_label(
                        current_final_count=current_final_count,
                        potential_additional=potential_additional,
                        component_count=component_count,
                    )
                    if potential_additional > 0
                    else "none"
                ),
                "interpretation": (
                    "candidate_for_stricter_ME_merge"
                    if potential_additional > 0
                    else "remains_one_component_at_0p10"
                ),
            }
        )

    candidate_rows = [item for item in module_rows if item["candidate"]]
    potential_additional = int(
        sum(item["potential_additional_modules"] for item in candidate_rows)
    )
    priority_order = {"high": 3, "medium": 2, "low": 1}
    network_priority = max(
        (item["priority"] for item in candidate_rows),
        key=lambda value: priority_order[value],
        default="none",
    )
    candidate_final_count = current_final_count + potential_additional
    figure_path = output_root / "figures" / safe_filename(key)
    draw_diagnostic_figure(
        linkage_matrix=linkage_matrix,
        labels=labels,
        figure_path=figure_path,
        title=(
            f"{clean_text(row['cell_type'])} | {clean_text(row['specimen'])} | "
            f"{clean_text(row['threshold'])} | {clean_text(row['strategy'])}"
        ),
        current_final_count=current_final_count,
        candidate_module_count=candidate_final_count,
        potential_additional=potential_additional,
    )

    actual_root = root / "05_find_modules" / ACTUAL_CONDITION
    actual_summary_path = actual_root / "analysis_summary.json"
    actual: dict[str, Any] = {}
    if actual_summary_path.exists():
        candidate_actual = json.loads(actual_summary_path.read_text(encoding="utf-8"))
        if candidate_actual.get("status") == "complete_findModules_merge_recalculation":
            actual = candidate_actual

    network_row = {
        "network_key": key,
        "threshold": clean_text(row["threshold"]),
        "strategy": clean_text(row["strategy"]),
        "mode": clean_text(row["mode"]),
        "specimen": clean_text(row["specimen"]),
        "cell_type": clean_text(row["cell_type"]),
        "samples": int(analysis["samples"]),
        "genes": int(analysis["genes"]),
        "selected_power": int(analysis["selected_power"]),
        "deepSplit": int(analysis["deepSplit"]),
        "minModuleSize": int(analysis["minModuleSize"]),
        "current_MEDissThres": float(analysis["MEDissThres"]),
        "candidate_MEDissThres": CANDIDATE_THRESHOLD,
        "dynamic_module_count": int(analysis["dynamic_module_count"]),
        "current_final_module_count": current_final_count,
        "diagnostic_tree_components_at_0p20": int(len(set(clusters_020))),
        "diagnostic_tree_components_at_0p10": int(len(set(clusters_010))),
        "candidate_current_module_count": int(len(candidate_rows)),
        "potential_additional_modules": potential_additional,
        "diagnostic_candidate_module_count": candidate_final_count,
        "collapsed_at_current_setting": bool(
            analysis["collapsed_to_one_proper_module"]
        ),
        "priority": network_priority,
        "recommended_first_test": (
            "rerun_MEDissThres_0.10" if candidate_rows else "keep_MEDissThres_0.20"
        ),
        "deepSplit_assessment": "not_inferred_from_ME_dendrogram",
        "actual_rerun_status": actual.get("status"),
        "actual_final_module_count": actual.get("actual_final_module_count"),
        "actual_additional_modules": actual.get("actual_additional_modules"),
        "prediction_minus_actual": actual.get("prediction_minus_actual"),
        "actual_collapsed": actual.get("collapsed_to_one_proper_module"),
        "diagnostic_figure": relative_link(figure_path, output_root),
        "source_me_dendrogram": relative_link(
            module_root / "figures/module_eigengene_dendrogram.png", output_root
        ),
        "source_analysis_summary": relative_link(analysis_path, output_root),
        "actual_gene_dendrogram": (
            relative_link(
                actual_root / "figures/gene_dendrogram_modules.png", output_root
            )
            if actual
            else None
        ),
        "actual_module_sizes": (
            relative_link(actual_root / "figures/module_sizes.png", output_root)
            if actual
            else None
        ),
        "actual_analysis_summary": (
            relative_link(actual_summary_path, output_root) if actual else None
        ),
    }
    actual_inputs = [actual_summary_path] if actual else []
    return network_row, module_rows, [*required, *actual_inputs]


def html_page(network_rows: list[dict[str, Any]], module_rows: list[dict[str, Any]]) -> str:
    records = []
    candidates_by_network: dict[str, list[dict[str, Any]]] = {}
    for item in module_rows:
        if item["candidate"]:
            candidates_by_network.setdefault(item["network_key"], []).append(item)
    for row in network_rows:
        record = dict(row)
        record["modules"] = candidates_by_network.get(row["network_key"], [])
        records.append(record)
    data = json.dumps(records, ensure_ascii=False).replace("</", "<\\/")
    return f"""<!doctype html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ME merge candidates</title>
<style>
:root{{--surface:#fff;--text:#262626;--muted:#666;--line:#ddd;--subtle:#eee;--accent:#007c83;--warn:#b14e23;--font:Arial,Helvetica,"Malgun Gothic",sans-serif}}
*{{box-sizing:border-box}} body{{margin:0;background:var(--surface);color:var(--text);font:13px/1.55 var(--font)}}
header{{max-width:1480px;margin:32px auto 0;padding:0 28px 24px;border-bottom:1px solid var(--subtle)}}
h1{{font-size:24px;font-weight:400;margin:0 0 8px}} .subtitle{{color:var(--muted);max-width:1100px}}
.method{{margin-top:12px;color:#444}} code{{font-family:Consolas,monospace;font-size:12px}}
.controls{{position:sticky;top:0;z-index:4;background:white;border-bottom:1px solid var(--line);padding:12px max(28px,calc((100vw - 1480px)/2 + 28px));display:flex;gap:8px;flex-wrap:wrap;align-items:center}}
select,input,button{{font:inherit;border:1px solid #bbb;background:white;padding:7px 9px;color:var(--text)}} button{{cursor:pointer}} button:hover{{background:#f5f5f5}} .count{{margin-left:auto;font-variant-numeric:tabular-nums}}
main{{max-width:1480px;margin:0 auto;padding:24px 28px 48px}} article{{padding:0 0 36px;margin:0 0 36px;border-bottom:1px solid var(--line)}}
h2{{font-size:16px;font-weight:500;margin:0 0 6px}} .facts{{display:flex;gap:6px 16px;flex-wrap:wrap;color:var(--muted);margin-bottom:12px}}
.priority{{font-weight:700;text-transform:uppercase}} .high{{color:#a33a19}} .medium{{color:#8b6500}} .low{{color:#386b73}}
.figure-grid{{display:grid;grid-template-columns:1fr 1fr;gap:16px}} .figure-wrap{{min-width:0}} .figure-title{{font-size:12px;color:var(--muted);margin:0 0 6px}} .figure{{display:block;border-top:1px solid var(--subtle);border-bottom:1px solid var(--subtle);overflow:auto;background:white}} .figure img{{display:block;width:100%;height:auto;min-height:240px;object-fit:contain}}
table{{border-collapse:collapse;width:100%;margin-top:14px;font-size:12px}} caption{{text-align:left;font-size:14px;font-weight:500;padding:0 0 8px}} th{{font-weight:400;color:var(--muted);border-bottom:1px solid #bbb;text-align:left;padding:7px 10px 7px 0}} td{{border-bottom:1px solid var(--subtle);padding:7px 10px 7px 0;vertical-align:top}} .number{{text-align:right;font-variant-numeric:tabular-nums}} a{{color:#444;text-underline-offset:3px}} .swatch{{display:inline-block;width:10px;height:10px;border:1px solid #999;margin-right:6px;vertical-align:-1px}}
.decision{{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin-top:12px}} .decision label{{white-space:nowrap}} .decision input[type=radio]{{padding:0}} .notes{{min-width:260px;flex:1}}
.empty{{color:var(--muted);padding:20px 0}} @media(max-width:900px){{.figure-grid{{grid-template-columns:1fr}}}} @media(max-width:700px){{header,main{{padding-left:16px;padding-right:16px}}.count{{width:100%;margin-left:0}}}}
@media print{{.controls,.decision{{display:none}}}}
</style>
</head>
<body>
<header>
<h1>ME dendrogram merge candidates</h1>
<div class="subtitle">현재 <code>MEDissThres=0.20</code>, 단순 <code>0.10</code> 절단 예측, PyWGCNA 공식 반복 병합을 <code>0.10</code>으로 실제 재계산한 결과를 함께 비교합니다.</div>
<div class="method"><strong>해석:</strong> candidate는 반복 ME 재계산 전 진단값이고, actual final이 공식 <code>mergeCloseModules(iterate=True)</code> 재계산 결과입니다. <code>deepSplit</code>은 gene dendrogram 단계의 network-wide 설정이므로 이 ME dendrogram에서는 판정하지 않았습니다.</div>
</header>
<div class="controls">
<select id="priority"><option value="">All priorities</option><option value="high">High</option><option value="medium">Medium</option><option value="low">Low</option><option value="none">No candidate</option></select>
<select id="threshold"><option value="">All thresholds</option></select>
<select id="strategy"><option value="">All strategies</option></select>
<select id="mode"><option value="">All modes</option></select>
<select id="specimen"><option value="">All specimens</option></select>
<select id="cell"><option value="">All cell types</option></select>
<input id="search" placeholder="Search cell type">
<button id="export">Export review CSV</button>
<span class="count" id="count"></span>
</div>
<main id="results"></main>
<script>
const records={data};
const storageKey='me-merge-candidate-review-v1';
let reviews=JSON.parse(localStorage.getItem(storageKey)||'{{}}');
const esc=v=>String(v??'').replace(/[&<>"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}}[c]));
const fields={{threshold:'threshold',strategy:'strategy',mode:'mode',specimen:'specimen',cell:'cell_type'}};
for(const [id,key] of Object.entries(fields)){{const node=document.getElementById(id);[...new Set(records.map(r=>r[key]))].sort().forEach(v=>{{const o=document.createElement('option');o.value=v;o.textContent=v;node.appendChild(o)}});node.addEventListener('change',render)}}
document.getElementById('priority').addEventListener('change',render); document.getElementById('search').addEventListener('input',render);
function moduleTable(row){{if(!row.modules.length)return '<div class="empty">0.10 절단선에서 추가 분리 후보가 없습니다.</div>';return `<table><caption>Current final modules that separate at 0.10</caption><thead><tr><th>Current module</th><th class="number">Features</th><th class="number">Dynamic modules</th><th class="number">0.10 components</th><th class="number">Potential +</th><th>Dynamic module sizes</th><th class="number">Max ME distance</th></tr></thead><tbody>${{row.modules.map(m=>`<tr><td><span class="swatch" style="background:${{esc(m.current_final_module)}}"></span>${{esc(m.current_final_module)}}</td><td class="number">${{m.feature_count}}</td><td class="number">${{m.dynamic_module_count}}</td><td class="number">${{m.diagnostic_components_at_0p10}}</td><td class="number">${{m.potential_additional_modules}}</td><td>${{esc(m.dynamic_module_sizes)}}</td><td class="number">${{m.pairwise_me_dissimilarity_max==null?'—':Number(m.pairwise_me_dissimilarity_max).toFixed(3)}}</td></tr>`).join('')}}</tbody></table>`}}
function decisionUI(row){{const saved=reviews[row.network_key]||{{decision:'hold',notes:''}};return `<div class="decision" data-key="${{esc(row.network_key)}}"><strong>Review</strong>${{['use_0.10','hold','keep_0.20'].map(v=>`<label><input type="radio" name="d-${{esc(row.network_key)}}" value="${{v}}" ${{saved.decision===v?'checked':''}}> ${{v==='use_0.10'?'Use 0.10':v==='hold'?'Hold':'Keep 0.20'}}</label>`).join('')}}<input class="notes" value="${{esc(saved.notes)}}" placeholder="Notes"></div>`}}
function render(){{const q=document.getElementById('search').value.trim().toLowerCase();const filtered=records.filter(r=>{{if(document.getElementById('priority').value&&r.priority!==document.getElementById('priority').value)return false;for(const [id,key] of Object.entries(fields)){{const wanted=document.getElementById(id).value;if(wanted&&r[key]!==wanted)return false}}return !q||`${{r.cell_type}} ${{r.network_key}}`.toLowerCase().includes(q)}});document.getElementById('count').textContent=`${{filtered.length}} / ${{records.length}} networks`;document.getElementById('results').innerHTML=filtered.map(row=>{{const actual=row.actual_final_module_count==null?'not run':row.actual_final_module_count;const actualLinks=row.actual_analysis_summary?`<a href="${{esc(row.actual_analysis_summary)}}" target="_blank">actual summary</a><a href="${{esc(row.actual_module_sizes)}}" target="_blank">actual sizes</a>`:'';const actualFigure=row.actual_gene_dendrogram?`<div class="figure-wrap"><div class="figure-title">Actual modules after iterative merge at 0.10</div><a class="figure" href="${{esc(row.actual_gene_dendrogram)}}" target="_blank"><img loading="lazy" src="${{esc(row.actual_gene_dendrogram)}}" alt="Actual gene dendrogram modules"></a></div>`:'';return `<article><h2>${{esc(row.cell_type)}} · ${{esc(row.specimen)}} <span class="priority ${{row.priority}}">${{esc(row.priority)}}</span></h2><div class="facts"><span>${{esc(row.threshold)}}</span><span>${{esc(row.strategy)}}</span><span>${{esc(row.mode)}}</span><span>n=${{row.samples}}</span><span>genes=${{row.genes}}</span><span>dynamic=${{row.dynamic_module_count}}</span><span>current final=${{row.current_final_module_count}}</span><span>candidate=${{row.diagnostic_candidate_module_count}}</span><strong>actual final=${{actual}}</strong><span>actual +${{row.actual_additional_modules??'—'}}</span><span>prediction - actual=${{row.prediction_minus_actual??'—'}}</span><a href="${{esc(row.source_me_dendrogram)}}" target="_blank">current figure</a><a href="${{esc(row.source_analysis_summary)}}" target="_blank">current summary</a>${{actualLinks}}</div><div class="figure-grid"><div class="figure-wrap"><div class="figure-title">Diagnostic ME tree: current 0.20 and candidate 0.10 lines</div><a class="figure" href="${{esc(row.diagnostic_figure)}}" target="_blank"><img loading="lazy" src="${{esc(row.diagnostic_figure)}}" alt="ME merge diagnostic"></a></div>${{actualFigure}}</div>${{moduleTable(row)}}${{decisionUI(row)}}</article>`}}).join('');document.querySelectorAll('.decision').forEach(box=>{{const key=box.dataset.key;box.querySelectorAll('input[type=radio]').forEach(node=>node.addEventListener('change',()=>saveReview(key,box)));box.querySelector('.notes').addEventListener('input',()=>saveReview(key,box))}})}}
function saveReview(key,box){{reviews[key]={{decision:box.querySelector('input[type=radio]:checked').value,notes:box.querySelector('.notes').value,updated_at:new Date().toISOString()}};localStorage.setItem(storageKey,JSON.stringify(reviews))}}
document.getElementById('export').addEventListener('click',()=>{{const header=['network_key','decision','notes','updated_at'];const quote=v=>'"'+String(v??'').replaceAll('"','""')+'"';const rows=records.map(r=>{{const d=reviews[r.network_key]||{{decision:'hold',notes:'',updated_at:''}};return [r.network_key,d.decision,d.notes,d.updated_at].map(quote).join(',')}});const blob=new Blob([[header.join(','),...rows].join('\\n')],{{type:'text/csv;charset=utf-8'}});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='me_merge_candidate_review.csv';a.click();URL.revokeObjectURL(a.href)}});
render();
</script>
</body>
</html>
"""


def run(run_root: Path) -> None:
    stage_root = run_root / "experiments/04_variance_filtered/transcriptomics"
    output_root = stage_root / OUTPUT_DIRECTORY
    output_root.mkdir(parents=True, exist_ok=True)
    registry_path = stage_root / "find_modules_summary.csv"
    registry = pd.read_csv(registry_path)
    registry = registry.loc[registry["status"] == "complete_findModules"].copy()
    registry = registry.sort_values(
        ["strategy", "threshold", "mode", "specimen", "cell_type"], kind="stable"
    )

    network_rows: list[dict[str, Any]] = []
    module_rows: list[dict[str, Any]] = []
    input_paths: list[Path] = [registry_path]
    errors: list[dict[str, str]] = []
    for position, (_, row) in enumerate(registry.iterrows(), start=1):
        print(f"[{position}/{len(registry)}] {network_key(row)}", flush=True)
        try:
            network_row, current_modules, current_inputs = analyze_network(
                stage_root=stage_root,
                output_root=output_root,
                row=row,
            )
            network_rows.append(network_row)
            module_rows.extend(current_modules)
            input_paths.extend(current_inputs)
        except Exception as error:
            errors.append({"network_key": network_key(row), "error": repr(error)})
            print(f"  ERROR: {error!r}", flush=True)

    network_table = pd.DataFrame(network_rows)
    module_table = pd.DataFrame(module_rows)
    priority_rank = {"high": 0, "medium": 1, "low": 2, "none": 3}
    if not network_table.empty:
        network_table["_priority"] = network_table["priority"].map(priority_rank)
        network_table = network_table.sort_values(
            ["_priority", "potential_additional_modules", "network_key"],
            ascending=[True, False, True],
            kind="stable",
        ).drop(columns="_priority")
    if not module_table.empty:
        module_table["_priority"] = module_table["priority"].map(priority_rank)
        module_table = module_table.sort_values(
            ["candidate", "_priority", "potential_additional_modules", "network_key"],
            ascending=[False, True, False, True],
            kind="stable",
        ).drop(columns="_priority")

    network_csv = output_root / "network_candidates.csv"
    module_csv = output_root / "module_candidates.csv"
    candidate_csv = output_root / "candidate_modules_only.csv"
    dashboard_path = output_root / "index.html"
    error_path = output_root / "errors.json"
    summary_path = output_root / "analysis_summary.json"
    manifest_path = output_root / "run_manifest.json"
    network_table.to_csv(network_csv, index=False)
    module_table.to_csv(module_csv, index=False)
    module_table.loc[module_table["candidate"]].to_csv(candidate_csv, index=False)
    dashboard_path.write_text(
        html_page(network_table.to_dict("records"), module_table.to_dict("records")),
        encoding="utf-8",
    )
    write_json(error_path, {"errors": errors})

    candidate_networks = network_table.loc[
        network_table["potential_additional_modules"] > 0
    ]
    candidate_modules = module_table.loc[module_table["candidate"]]
    summary = {
        "status": "complete" if not errors else "complete_with_errors",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_type": "diagnostic_only_no_findModules_rerun",
        "baseline_MEDissThres": BASELINE_THRESHOLD,
        "candidate_MEDissThres": CANDIDATE_THRESHOLD,
        "networks_analyzed": int(len(network_table)),
        "networks_with_candidate": int(len(candidate_networks)),
        "current_final_modules_with_candidate": int(len(candidate_modules)),
        "high_priority_networks": int((network_table["priority"] == "high").sum()),
        "medium_priority_networks": int(
            (network_table["priority"] == "medium").sum()
        ),
        "low_priority_networks": int((network_table["priority"] == "low").sum()),
        "collapsed_networks_with_candidate": int(
            (
                candidate_networks["collapsed_at_current_setting"].astype(bool)
            ).sum()
        ),
        "total_potential_additional_modules": int(
            candidate_modules["potential_additional_modules"].sum()
        ),
        "errors": int(len(errors)),
        "deepSplit_note": (
            "deepSplit acts on the gene dendrogram before ME merging and cannot be "
            "selected from this ME dendrogram diagnostic"
        ),
        "interpretation_limit": (
            "The 0.10 components are pre-rerun candidates, not exact final modules, "
            "because PyWGCNA iteratively recomputes eigengenes during merging"
        ),
    }
    write_json(summary_path, summary)
    outputs = [
        network_csv,
        module_csv,
        candidate_csv,
        dashboard_path,
        error_path,
        summary_path,
        *sorted((output_root / "figures").glob("*.png")),
    ]
    write_run_manifest(
        manifest_path,
        pipeline_name="transcriptomics_me_merge_candidate_diagnostic",
        run_root=run_root,
        parameters={
            "baseline_MEDissThres": BASELINE_THRESHOLD,
            "candidate_MEDissThres": CANDIDATE_THRESHOLD,
            "linkage": "saved average-linkage tree",
            "distance": "1 - Pearson correlation of dynamic module eigengenes",
            "findModules_rerun": False,
        },
        generation_rules={
            "candidate": (
                "a current final module whose proper dynamic modules occupy more "
                "than one saved-tree component at MEDissThres=0.10"
            ),
            "raw_dynamic_label_0": "excluded as unassigned",
            "deepSplit": "not inferred because it acts on the gene dendrogram",
            "exactness": "diagnostic only; exact final assignments require rerun",
        },
        input_paths=input_paths,
        output_paths=outputs,
        code_paths=[Path(__file__), Path(__file__).parent / "pipeline_provenance.py"],
        packages=("numpy", "pandas", "matplotlib", "scipy"),
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    run(arguments.run_root.resolve())
