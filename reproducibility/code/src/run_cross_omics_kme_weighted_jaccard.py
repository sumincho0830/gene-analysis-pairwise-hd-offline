"""Compare every transcriptomics module with specimen-matched proteomics modules.

The descriptive similarity is the established project definition

    sum_g min(|kME_RNA(g)|, |kME_protein(g)|)
    ------------------------------------------------
    sum_g max(|kME_RNA(g)|, |kME_protein(g)|)

where each weight is non-zero only for genes assigned to that module.  The
universe is the intersection of genes measurable in both assays for the given
transcriptomics experiment.  Raw-label-0 unassigned groups are excluded.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import seaborn as sns

from pipeline_provenance import write_run_manifest


TRANSCRIPT_CONDITION = "01_default_ds2_mms50_med0p20"
PROTEOMICS_REFERENCE_CONDITION = "04_ds3_mms30_med0p10"
AMBIGUOUS_SYMBOL = re.compile(r"[|;,/]")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def normalize_gene(value: Any) -> str | None:
    if pd.isna(value):
        return None
    gene = str(value).strip()
    if not gene or gene.lower() in {"nan", "none", "na"}:
        return None
    if AMBIGUOUS_SYMBOL.search(gene):
        return None
    return gene.upper()


def weighted_jaccard(
    left: dict[str, float], right: dict[str, float]
) -> tuple[float, float, int, float, float, list[dict[str, Any]]]:
    universe = sorted(set(left) | set(right))
    if not universe:
        return np.nan, np.nan, 0, 0.0, 0.0, []
    left_values = np.asarray([left.get(gene, 0.0) for gene in universe], dtype=float)
    right_values = np.asarray([right.get(gene, 0.0) for gene in universe], dtype=float)
    minimum = np.minimum(left_values, right_values)
    maximum = np.maximum(left_values, right_values)
    numerator = float(minimum.sum())
    denominator = float(maximum.sum())
    shared = set(left) & set(right)
    plain = len(shared) / len(universe)
    contributions = [
        {
            "gene": gene,
            "transcript_abs_kME": float(left.get(gene, 0.0)),
            "proteomics_abs_kME": float(right.get(gene, 0.0)),
            "weighted_intersection": float(minimum[index]),
            "weighted_union": float(maximum[index]),
        }
        for index, gene in enumerate(universe)
        if gene in shared
    ]
    return (
        numerator / denominator if denominator else np.nan,
        float(plain),
        len(shared),
        numerator,
        denominator,
        contributions,
    )


def collapse_proteomics_to_gene(
    table: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str], dict[str, Any]]:
    """Select one informative proper-module probe for each simple gene symbol.

    This follows the purpose of WGCNA ``collapseRowsUsingKME``: duplicate
    probes are represented by the probe with the largest signed kME to any
    proper module.  Ties are resolved by PROBEID for determinism.
    """

    data = table.copy()
    data["gene_symbol"] = data["SYMBOL"].map(normalize_gene)
    valid_simple = data["gene_symbol"].notna()
    universe = set(data.loc[valid_simple, "gene_symbol"].astype(str))
    ambiguous = data.loc[data["SYMBOL"].astype(str).str.contains(AMBIGUOUS_SYMBOL, na=False)]
    proper = data.loc[
        valid_simple
        & ~data["is_raw_unassigned"].astype(bool)
        & data["official_assigned_abs_kME"].notna()
    ].copy()
    kme_columns = [column for column in proper.columns if column.startswith("kME")]
    if not kme_columns:
        raise ValueError("Proteomics table has no proper-module kME columns")
    proper["selection_max_signed_kME_any_module"] = proper[kme_columns].max(axis=1)
    probe_counts = proper.groupby("gene_symbol").size().rename("proper_probe_count")
    module_counts = (
        proper.groupby("gene_symbol")["module_color"]
        .nunique()
        .rename("proper_probe_module_count")
    )
    proper = proper.sort_values(
        ["gene_symbol", "selection_max_signed_kME_any_module", "PROBEID"],
        ascending=[True, False, True],
        kind="stable",
    )
    selected = proper.groupby("gene_symbol", sort=True, as_index=False).head(1).copy()
    selected = selected.join(probe_counts, on="gene_symbol").join(
        module_counts, on="gene_symbol"
    )
    selected = selected[
        [
            "gene_symbol",
            "PROBEID",
            "SYMBOL",
            "module_color",
            "official_assigned_signed_kME",
            "official_assigned_abs_kME",
            "selection_max_signed_kME_any_module",
            "proper_probe_count",
            "proper_probe_module_count",
        ]
    ].sort_values("gene_symbol", kind="stable")
    audit = {
        "proteomics_features": int(len(data)),
        "simple_symbol_features": int(valid_simple.sum()),
        "ambiguous_symbol_features_excluded": int(len(ambiguous)),
        "unique_simple_symbol_universe": int(len(universe)),
        "proper_module_candidate_probes": int(len(proper)),
        "selected_representative_genes": int(len(selected)),
        "genes_with_multiple_proper_probes": int((probe_counts > 1).sum()),
        "genes_with_probes_in_multiple_modules": int((module_counts > 1).sum()),
        "selection_rule": "highest signed kME across all proper modules; deterministic PROBEID tie-break",
    }
    return selected, universe, audit


def module_weights(
    table: pd.DataFrame,
    gene_column: str,
    weight_column: str,
    common_universe: set[str],
    module_names: list[str] | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, int]]:
    proper_mask = (
        ~table["is_raw_unassigned"].astype(bool)
        & table[weight_column].notna()
    )
    if module_names is None:
        module_names = list(
            dict.fromkeys(table.loc[proper_mask, "module_color"].astype(str).tolist())
        )
    proper = table.loc[proper_mask & table[gene_column].isin(common_universe)].copy()
    output: dict[str, dict[str, float]] = {str(module): {} for module in module_names}
    counts: dict[str, int] = {str(module): 0 for module in module_names}
    for module, frame in proper.groupby("module_color", sort=False):
        collapsed = frame.groupby(gene_column, sort=False)[weight_column].max()
        output[str(module)] = collapsed.astype(float).to_dict()
        counts[str(module)] = int(len(collapsed))
    return output, counts


def compare_modules(
    transcript_modules: dict[str, dict[str, float]],
    protein_modules: dict[str, dict[str, float]],
    transcript_counts: dict[str, int],
    protein_counts: dict[str, int],
    common_universe_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    contribution_cache: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for transcript_module, left in transcript_modules.items():
        for protein_module, right in protein_modules.items():
            weighted, plain, shared, numerator, denominator, contributions = weighted_jaccard(
                left, right
            )
            contribution_cache[(transcript_module, protein_module)] = contributions
            rows.append(
                {
                    "transcript_module": transcript_module,
                    "proteomics_module": protein_module,
                    "common_gene_universe": common_universe_size,
                    "transcript_module_genes_in_common": transcript_counts[transcript_module],
                    "proteomics_module_genes_in_common": protein_counts[protein_module],
                    "shared_genes": shared,
                    "jaccard": plain,
                    "weighted_intersection": numerator,
                    "weighted_union": denominator,
                    "kme_weighted_jaccard": weighted,
                }
            )
    result = pd.DataFrame(rows)
    if result.empty:
        return result, pd.DataFrame()
    result["transcript_best"] = (
        result["transcript_module_genes_in_common"].gt(0)
        & result["kme_weighted_jaccard"].eq(
            result.groupby("transcript_module")["kme_weighted_jaccard"].transform("max")
        )
    )
    result["proteomics_best"] = (
        result["proteomics_module_genes_in_common"].gt(0)
        & result["kme_weighted_jaccard"].eq(
            result.groupby("proteomics_module")["kme_weighted_jaccard"].transform("max")
        )
    )
    result["reciprocal_best"] = result["transcript_best"] & result["proteomics_best"]
    result["transcript_rank"] = result.groupby("transcript_module")[
        "kme_weighted_jaccard"
    ].rank(method="min", ascending=False)
    result["proteomics_rank"] = result.groupby("proteomics_module")[
        "kme_weighted_jaccard"
    ].rank(method="min", ascending=False)

    contribution_rows: list[dict[str, Any]] = []
    for (transcript_module, proteomics_module), items in contribution_cache.items():
        for item in items:
            contribution_rows.append(
                {
                    "transcript_module": transcript_module,
                    "proteomics_module": proteomics_module,
                    **item,
                }
            )
    return result, pd.DataFrame(
        contribution_rows,
        columns=[
            "transcript_module",
            "proteomics_module",
            "gene",
            "transcript_abs_kME",
            "proteomics_abs_kME",
            "weighted_intersection",
            "weighted_union",
        ],
    )


def plot_similarity(
    table: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    subtitle: str,
    global_vmax: float,
) -> None:
    transcript_order = (
        table[["transcript_module", "transcript_module_genes_in_common"]]
        .drop_duplicates()
        .sort_values(
            ["transcript_module_genes_in_common", "transcript_module"],
            ascending=[False, True],
            kind="stable",
        )["transcript_module"]
        .tolist()
    )
    protein_order = (
        table[["proteomics_module", "proteomics_module_genes_in_common"]]
        .drop_duplicates()
        .sort_values(
            ["proteomics_module_genes_in_common", "proteomics_module"],
            ascending=[False, True],
            kind="stable",
        )["proteomics_module"]
        .tolist()
    )
    matrix = table.pivot(
        index="transcript_module",
        columns="proteomics_module",
        values="kme_weighted_jaccard",
    ).reindex(index=transcript_order, columns=protein_order)
    transcript_counts = (
        table.set_index("transcript_module")["transcript_module_genes_in_common"]
        .to_dict()
    )
    protein_counts = (
        table.set_index("proteomics_module")["proteomics_module_genes_in_common"]
        .to_dict()
    )
    width = max(11.0, len(protein_order) * 0.58 + 4.5)
    height = max(6.5, len(transcript_order) * 0.38 + 3.5)
    figure, axis = plt.subplots(figsize=(width, height), facecolor="white")
    sns.heatmap(
        matrix,
        cmap="mako",
        vmin=0,
        vmax=global_vmax,
        linewidths=0.25,
        linecolor="#f2f4f7",
        cbar_kws={"label": "|kME|-weighted Jaccard (global scale)"},
        ax=axis,
    )
    axis.set_xticklabels(
        [f"{module}\n(n={protein_counts[module]:,})" for module in protein_order],
        rotation=55,
        ha="right",
        fontsize=8,
    )
    axis.set_yticklabels(
        [f"{module} (n={transcript_counts[module]:,})" for module in transcript_order],
        rotation=0,
        fontsize=8,
    )
    positions_y = {module: index for index, module in enumerate(transcript_order)}
    positions_x = {module: index for index, module in enumerate(protein_order)}
    row_best = table.loc[table["transcript_best"]]
    for match in row_best.itertuples(index=False):
        x = positions_x[match.proteomics_module]
        y = positions_y[match.transcript_module]
        marker = "★" if match.reciprocal_best else "•"
        value = float(match.kme_weighted_jaccard)
        text_color = "white" if value >= global_vmax * 0.48 else "black"
        axis.text(
            x + 0.5,
            y + 0.5,
            f"{value:.2f}{marker}",
            ha="center",
            va="center",
            fontsize=7,
            fontweight="bold",
            color=text_color,
        )
        if match.reciprocal_best:
            axis.add_patch(
                Rectangle((x, y), 1, 1, fill=False, edgecolor="#f59e0b", linewidth=2.0)
            )
    axis.set_xlabel("Proteomics proper modules (representative gene count in common universe)")
    axis.set_ylabel("Transcriptomics proper modules (gene count in common universe)")
    figure.suptitle(title, fontsize=14, fontweight="bold", y=0.995)
    axis.set_title(
        subtitle + "\n• RNA-module best; ★ reciprocal best (orange border)",
        fontsize=9,
        color="#475569",
        pad=12,
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def network_root(stage_root: Path, row: Any) -> Path:
    root = stage_root / row.threshold / row.strategy / row.mode / row.specimen
    if row.mode == "cell_separate":
        root = root / row.cell_type
    return root


def output_root(cross_root: Path, row: Any) -> Path:
    root = (
        cross_root
        / "transcriptomics"
        / row.threshold
        / row.strategy
        / row.mode
        / row.specimen
    )
    if row.mode == "cell_separate":
        root = root / row.cell_type
    return root


def build_gallery(summary: pd.DataFrame, cross_root: Path) -> None:
    records = []
    for row in summary.itertuples(index=False):
        records.append(
            {
                "threshold": row.threshold,
                "strategy": row.strategy,
                "mode": row.mode,
                "specimen": row.specimen,
                "cell_type": row.cell_type,
                "network_id": row.network_id,
                "common_genes": int(row.common_gene_universe),
                "rna_modules": int(row.transcript_module_count),
                "protein_modules": int(row.proteomics_module_count),
                "max_similarity": float(row.maximum_weighted_jaccard),
                "reciprocal": int(row.reciprocal_best_pairs),
                "image": Path(row.figure).as_posix(),
                "table": Path(row.pair_table).as_posix(),
            }
        )
    payload = json.dumps(records, ensure_ascii=False, separators=(",", ":"))
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>RNA-protein module similarity</title>
<style>body{{font-family:Arial,sans-serif;margin:0;background:#f4f6f8;color:#172033}}header{{padding:18px 24px;background:white;border-bottom:1px solid #d8dee8;position:sticky;top:0;z-index:3}}h1{{margin:0 0 10px;font-size:22px}}.filters{{display:flex;gap:8px;flex-wrap:wrap}}select,input{{padding:7px;border:1px solid #bdc7d6;border-radius:5px;background:white}}#count{{margin-left:auto;color:#526178}}main{{display:grid;grid-template-columns:repeat(auto-fill,minmax(500px,1fr));gap:14px;padding:16px}}.card{{background:white;border:1px solid #d8dee8;border-radius:8px;padding:12px}}.title{{font-weight:700}}.meta{{font-size:12px;color:#607089;margin:6px 0}}img{{width:100%;height:390px;object-fit:contain;background:white}}a{{color:#087f8c}}@media(max-width:650px){{main{{grid-template-columns:1fr}}}}</style></head><body>
<header><h1>RNA-protein kME-weighted Jaccard</h1><div class='filters'><select id='threshold'><option value=''>All thresholds</option></select><select id='strategy'><option value=''>All strategies</option></select><select id='mode'><option value=''>All modes</option></select><select id='specimen'><option value=''>All specimens</option></select><select id='cell'><option value=''>All cell types</option></select><input id='search' placeholder='Search network'><span id='count'></span></div></header><main id='grid'></main>
<script>const data={payload};const $=x=>document.getElementById(x);function values(k){{return [...new Set(data.map(x=>x[k]))].sort()}}function fill(id,k){{for(const v of values(k)){{const o=document.createElement('option');o.value=v;o.textContent=v;$(id).appendChild(o)}}$(id).onchange=render}}fill('threshold','threshold');fill('strategy','strategy');fill('mode','mode');fill('specimen','specimen');fill('cell','cell_type');$('search').oninput=render;function render(){{const q=$('search').value.toLowerCase();const rows=data.filter(x=>(!$('threshold').value||x.threshold===$('threshold').value)&&(!$('strategy').value||x.strategy===$('strategy').value)&&(!$('mode').value||x.mode===$('mode').value)&&(!$('specimen').value||x.specimen===$('specimen').value)&&(!$('cell').value||x.cell_type===$('cell').value)&&(!q||x.network_id.toLowerCase().includes(q)));$('count').textContent=`${{rows.length}} / ${{data.length}}`;$('grid').innerHTML=rows.map(x=>`<article class='card'><div class='title'>${{x.cell_type}} · ${{x.specimen}}</div><div class='meta'>${{x.threshold}} · ${{x.strategy}} · ${{x.mode}} · common genes ${{x.common_genes.toLocaleString()}} · modules ${{x.rna_modules}}×${{x.protein_modules}} · max ${{x.max_similarity.toFixed(3)}} · reciprocal ${{x.reciprocal}}</div><a href='${{x.image}}' target='_blank'><img loading='lazy' src='${{x.image}}'></a><div><a href='${{x.table}}'>pair table CSV</a></div></article>`).join('')}}render();</script></body></html>"""
    (cross_root / "kme_weighted_jaccard_gallery.html").write_text(html, encoding="utf-8")


def run(run_root: Path, branch: str = "baseline") -> None:
    transcript_stage = run_root / "experiments/04_variance_filtered/transcriptomics"
    cross_root = run_root / "experiments/05_cross_omics_module_similarity"
    after_outlier_removal = branch == "after_outlier_removal"
    output_stage = cross_root / "after_outlier_removal" if after_outlier_removal else cross_root
    shared_root = cross_root / "00_shared/proteomics_gene_collapse"
    shared_root.mkdir(parents=True, exist_ok=True)

    protein_cache: dict[str, dict[str, Any]] = {}
    for specimen in ("CSF", "PB"):
        source = (
            run_root
            / "experiments/01_batch_corrected/proteomics/06_kme_module_diagnosis"
            / specimen
            / PROTEOMICS_REFERENCE_CONDITION
            / "feature_module_kme.csv.gz"
        )
        table = pd.read_csv(source)
        collapsed, universe, audit = collapse_proteomics_to_gene(table)
        collapsed_path = shared_root / f"{specimen}_representative_gene_modules.csv.gz"
        collapsed.to_csv(
            collapsed_path,
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        write_json(shared_root / f"{specimen}_collapse_audit.json", audit)
        protein_cache[specimen] = {
            "source": source,
            "collapsed": collapsed,
            "universe": universe,
            "audit": audit,
            "collapsed_path": collapsed_path,
            "proper_module_names": list(
                dict.fromkeys(
                    table.loc[
                        ~table["is_raw_unassigned"].astype(bool)
                        & table["official_assigned_abs_kME"].notna(),
                        "module_color",
                    ]
                    .astype(str)
                    .tolist()
                )
            ),
        }

    registry = pd.read_csv(
        transcript_stage
        / (
            "outlier_removal_find_modules_summary.csv"
            if after_outlier_removal
            else "find_modules_summary.csv"
        )
    )
    registry = registry.loc[registry["status"].eq("complete_findModules")].copy()
    completed: list[dict[str, Any]] = []
    all_values: list[float] = []
    print(f"Computing {len(registry)} transcriptomics-to-proteomics module comparisons", flush=True)
    for row in registry.itertuples(index=False):
        source_root = network_root(transcript_stage, row)
        if after_outlier_removal:
            source_root = source_root / "outlier_removal"
        transcript_path = (
            source_root
            / "06_kme_module_diagnosis"
            / TRANSCRIPT_CONDITION
            / "feature_module_kme.csv.gz"
        )
        transcript = pd.read_csv(transcript_path)
        transcript["gene_symbol"] = transcript["gene_id"].map(normalize_gene)
        transcript_universe = set(transcript["gene_symbol"].dropna().astype(str))
        protein = protein_cache[row.specimen]
        common = transcript_universe & protein["universe"]

        protein_collapsed = protein["collapsed"].copy()
        protein_collapsed["is_raw_unassigned"] = False
        transcript_module_names = list(
            dict.fromkeys(
                transcript.loc[
                    ~transcript["is_raw_unassigned"].astype(bool)
                    & transcript["official_assigned_abs_kME"].notna(),
                    "module_color",
                ]
                .astype(str)
                .tolist()
            )
        )
        transcript_modules, transcript_counts = module_weights(
            transcript,
            "gene_symbol",
            "official_assigned_abs_kME",
            common,
            transcript_module_names,
        )
        protein_modules, protein_counts = module_weights(
            protein_collapsed,
            "gene_symbol",
            "official_assigned_abs_kME",
            common,
            protein["proper_module_names"],
        )
        pairs, all_contributions = compare_modules(
            transcript_modules,
            protein_modules,
            transcript_counts,
            protein_counts,
            len(common),
        )
        out = output_root(output_stage, row)
        out.mkdir(parents=True, exist_ok=True)
        pair_path = out / "module_pair_similarity.csv"
        top_path = out / "transcript_module_best_matches.csv"
        all_contribution_path = out / "all_pair_gene_contributions.csv.gz"
        contribution_path = out / "reciprocal_best_gene_contributions.csv.gz"
        pairs.to_csv(pair_path, index=False)
        top = pairs.loc[pairs["transcript_best"]].sort_values(
            ["transcript_module", "kme_weighted_jaccard", "proteomics_module"],
            ascending=[True, False, True],
            kind="stable",
        )
        top.to_csv(top_path, index=False)
        all_contributions.to_csv(
            all_contribution_path,
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        reciprocal_pairs = pairs.loc[
            pairs["reciprocal_best"], ["transcript_module", "proteomics_module"]
        ]
        contributions = all_contributions.merge(
            reciprocal_pairs,
            on=["transcript_module", "proteomics_module"],
            how="inner",
        )
        contributions.to_csv(
            contribution_path,
            index=False,
            compression={"method": "gzip", "compresslevel": 6, "mtime": 0},
        )
        mapping_audit = {
            "transcript_feature_universe": len(transcript_universe),
            "proteomics_simple_symbol_universe": len(protein["universe"]),
            "common_gene_universe": len(common),
            "transcript_proper_modules": len(transcript_modules),
            "proteomics_proper_modules": len(protein_modules),
            "transcript_unassigned_features_excluded_from_modules": int(
                transcript["is_raw_unassigned"].astype(bool).sum()
            ),
            "proteomics_reference_condition": PROTEOMICS_REFERENCE_CONDITION,
            "gene_normalization": "trim and uppercase",
            "ambiguous_proteomics_symbols": "symbols containing | ; , or / excluded from gene mapping",
            "comparison_universe": "intersection of transcript feature genes and all simple proteomics symbols",
        }
        audit_path = out / "gene_mapping_audit.json"
        write_json(audit_path, mapping_audit)
        figure_path = out / "kme_weighted_jaccard_heatmap.png"
        all_values.extend(
            pairs["kme_weighted_jaccard"].dropna().astype(float).tolist()
        )
        completed.append(
            {
                "row": row,
                "source_root": source_root,
                "transcript_path": transcript_path,
                "protein_source": protein["source"],
                "protein_collapsed_path": protein["collapsed_path"],
                "out": out,
                "pairs": pairs,
                "pair_path": pair_path,
                "top_path": top_path,
                "all_contribution_path": all_contribution_path,
                "all_contribution_rows": int(len(all_contributions)),
                "contribution_path": contribution_path,
                "audit_path": audit_path,
                "figure_path": figure_path,
                "common": len(common),
            }
        )

    if not all_values:
        raise RuntimeError("No finite weighted Jaccard values were calculated")
    global_vmax = float(max(all_values))
    summary_rows: list[dict[str, Any]] = []
    for index, item in enumerate(completed, 1):
        row = item["row"]
        pairs = item["pairs"]
        figure_path = item["figure_path"]
        plot_similarity(
            pairs,
            figure_path,
            title="Transcriptomics–proteomics module similarity",
            subtitle=(
                f"{row.threshold} · {row.strategy} · {row.mode} · {row.specimen} · {row.cell_type}; "
                f"protein={PROTEOMICS_REFERENCE_CONDITION}; common genes={item['common']:,}"
            ),
            global_vmax=global_vmax,
        )
        best = pairs.loc[pairs["transcript_best"], "kme_weighted_jaccard"]
        summary = {
            "threshold": row.threshold,
            "strategy": row.strategy,
            "mode": row.mode,
            "specimen": row.specimen,
            "cell_type": row.cell_type,
            "network_id": row.network_id,
            "proteomics_condition": PROTEOMICS_REFERENCE_CONDITION,
            "analysis_branch": branch,
            "common_gene_universe": item["common"],
            "transcript_module_count": int(pairs["transcript_module"].nunique()),
            "proteomics_module_count": int(pairs["proteomics_module"].nunique()),
            "module_pairs": int(len(pairs)),
            "maximum_weighted_jaccard": float(pairs["kme_weighted_jaccard"].max()),
            "median_transcript_best_weighted_jaccard": float(best.median()),
            "reciprocal_best_pairs": int(pairs["reciprocal_best"].sum()),
            "global_figure_vmax": global_vmax,
            "figure": figure_path.relative_to(output_stage).as_posix(),
            "pair_table": item["pair_path"].relative_to(output_stage).as_posix(),
            "top_matches": item["top_path"].relative_to(output_stage).as_posix(),
            "all_pair_gene_contributions": item["all_contribution_path"].relative_to(output_stage).as_posix(),
            "all_pair_contribution_rows": item["all_contribution_rows"],
            "mapping_audit": item["audit_path"].relative_to(output_stage).as_posix(),
        }
        summary_rows.append(summary)
        summary_path = item["out"] / "analysis_summary.json"
        write_json(summary_path, summary)
        write_run_manifest(
            item["out"] / "run_manifest.json",
            pipeline_name="cross_omics_kme_weighted_jaccard",
            run_root=run_root,
            parameters={
                "proteomics_condition": PROTEOMICS_REFERENCE_CONDITION,
                "specimen_matching": row.specimen,
                "gene_symbol_normalization": "trim and uppercase",
                "proteomics_duplicate_probe_collapse": "highest signed kME across all proper modules",
                "weight": "absolute own-module kME",
                "global_figure_vmax": global_vmax,
            },
            generation_rules={
                "formula": "sum(min(abs own-module kME RNA, abs own-module kME protein)) / sum(max(weights))",
                "module_outside_weight": 0,
                "unassigned": "raw-label-0 groups excluded",
                "universe": "intersection of genes measurable in both assays for this transcript experiment",
                "interpretation": "descriptive fuzzy-membership overlap; not formal WGCNA module preservation or disease association",
            },
            input_paths=[
                item["transcript_path"],
                item["protein_source"],
                item["protein_collapsed_path"],
            ],
            output_paths=[
                item["pair_path"],
                item["top_path"],
                item["all_contribution_path"],
                item["contribution_path"],
                item["audit_path"],
                item["figure_path"],
                summary_path,
            ],
            code_paths=[Path(__file__), Path(__file__).parent / "pipeline_provenance.py"],
            packages=("numpy", "pandas", "matplotlib", "seaborn"),
        )
        if index % 20 == 0 or index == len(completed):
            print(f"Rendered {index}/{len(completed)} heatmaps", flush=True)

    summary_table = pd.DataFrame(summary_rows)
    summary_path = output_stage / "cross_omics_similarity_summary.csv"
    summary_table.to_csv(summary_path, index=False)
    write_json(
        output_stage / "method.json",
        {
            "method": "absolute-own-module-kME weighted Jaccard",
            "formula": "sum_g min(w_RNA,g, w_protein,g) / sum_g max(w_RNA,g, w_protein,g)",
            "weights": "absolute official assigned own-module kME; 0 outside module",
            "proteomics_reference_condition": PROTEOMICS_REFERENCE_CONDITION,
            "specimen_matching": "CSF RNA to CSF protein; PB RNA to PB protein",
            "universe": "per-experiment measurable gene intersection",
            "proteomics_probe_collapse": "WGCNA collapseRowsUsingKME principle; highest signed kME across proper modules",
            "global_color_scale_max": global_vmax,
            "interpretation_limit": "descriptive module membership overlap, not formal preservation or significance test",
            "analysis_branch": branch,
            "sample_exclusion": (
                "Ng_MOG1 and Ng_HD3; no finite sample-tree cut"
                if after_outlier_removal
                else "baseline stored branch"
            ),
        },
    )
    build_gallery(summary_table, output_stage)
    print(f"Completed {len(summary_table)} cross-omics experiments", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-root", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument(
        "--branch",
        choices=("baseline", "after_outlier_removal"),
        default="baseline",
    )
    args = parser.parse_args()
    run(args.run_root.resolve(), args.branch)
