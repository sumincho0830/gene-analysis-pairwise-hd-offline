(function () {
  const data = window.OUTLIER_CANDIDATE_DATA || {
    occurrences: [], samples: [], experiments: [], cell_type_experiment_sets: [], summary: {},
  };
  const groups = data.cell_type_experiment_sets || [];
  const storageKey = `outlier-review:${data.schema_version || '1'}:results_20261119`;
  let reviews = JSON.parse(localStorage.getItem(storageKey) || '{}');
  const state = {
    groupId: groups[0]?.id || '',
    query: '', decision: '', sort: 'weighted', selectedSample: null, recordId: null,
  };
  let searchTimer = null;

  const html = value => String(value ?? '').replace(/[&<>"']/g, char => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  }[char]));
  const finite = value => Number.isFinite(Number(value));
  const numberText = (value, digits = 2) => finite(value)
    ? Number(value).toLocaleString(undefined, {maximumFractionDigits: digits}) : 'NA';
  const percentileText = value => finite(value) ? `${numberText(value, 1)}%` : 'NA';
  const cellLabel = value => String(value || '').replaceAll('_', ' ');
  const reviewOf = sampleId => reviews[sampleId] || {
    decision: 'unreviewed', notes: '', updated_at: null,
  };
  const decisionLabel = value => ({
    keep: 'Keep', hold: 'Hold', exclude: 'Exclude', unreviewed: 'Unreviewed',
  }[value] || value);
  const hierarchyKeys = ['stage', 'specimen', 'threshold', 'strategy', 'phase'];
  const hierarchyMeta = {
    stage: ['Pipeline stage', 'Raw to variance filtered'],
    specimen: ['Specimen', 'ALL, CSF, or PB'],
    threshold: ['Threshold', '5 or 10 cells'],
    strategy: ['Feature strategy', 'CPM or HVG'],
    phase: ['Preprocess phase', 'Before or after'],
  };

  function currentGroup() {
    return groups.find(group => group.id === state.groupId) || groups[0] || null;
  }

  function persistReview(sampleId, patch) {
    reviews[sampleId] = {...reviewOf(sampleId), ...patch, updated_at: new Date().toISOString()};
    localStorage.setItem(storageKey, JSON.stringify(reviews));
    window.renderAll?.();
  }

  function visibleSamples(group = currentGroup()) {
    if (!group) return [];
    const query = state.query.trim().toLowerCase();
    const rows = group.candidate_samples.filter(row => {
      if (state.decision && reviewOf(row.source_sample_id).decision !== state.decision) return false;
      if (!query) return true;
      return `${row.source_sample_id} ${row.cohort} ${row.diagnosis} ${row.specimen}`
        .toLowerCase().includes(query);
    });
    const compare = {
      weighted: (a, b) => (b.pca_weighted_count ?? 0) - (a.pca_weighted_count ?? 0)
        || b.outlier_cell_count - a.outlier_cell_count
        || a.source_sample_id.localeCompare(b.source_sample_id),
      count: (a, b) => b.outlier_cell_count - a.outlier_cell_count
        || b.outlier_cell_fraction - a.outlier_cell_fraction
        || a.source_sample_id.localeCompare(b.source_sample_id),
      fraction: (a, b) => b.outlier_cell_fraction - a.outlier_cell_fraction
        || b.outlier_cell_count - a.outlier_cell_count
        || a.source_sample_id.localeCompare(b.source_sample_id),
      sample: (a, b) => a.source_sample_id.localeCompare(b.source_sample_id),
    }[state.sort];
    return [...rows].sort(compare);
  }

  function renderKpis() {
    const group = currentGroup();
    const rows = visibleSamples(group);
    const tested = rows.reduce((sum, row) => sum + row.tested_cell_count, 0);
    const flagged = rows.reduce((sum, row) => sum + row.outlier_cell_count, 0);
    const reviewed = rows.filter(row => reviewOf(row.source_sample_id).decision !== 'unreviewed').length;
    const values = [
      [groups.length, 'Experiment sets'],
      [group?.cell_types.length || 0, 'Cell types'],
      [rows.length, 'Candidate samples'],
      [flagged, 'Outlier cells'],
      [tested, 'Tested cells'],
      [reviewed, 'Reviewed'],
    ];
    document.getElementById('kpis').innerHTML = values.map(([value, label]) =>
      `<div class="kpi"><b>${numberText(value, 0)}</b><span>${label}</span></div>`).join('');
  }

  function hierarchyValueLabel(key, value) {
    const mappings = {
      stage: {
        '00_raw': 'Raw', '01_batch_corrected': 'Batch corrected',
        '02_sample_filtered': 'Sample filtered', '03_gene_filtered': 'Gene filtered',
        '04_variance_filtered': 'Variance filtered',
      },
      strategy: {
        '01_cpm_detection': 'CPM detection', '02_hvg_top9000': 'HVG top 9,000',
        '03_hvg_top5000': 'HVG top 5,000',
      },
      phase: {
        '02_data_distribution': 'Data distribution',
        '02_before_preprocess': 'Before preprocess',
        '04_after_preprocess': 'After preprocess',
      },
    };
    if (mappings[key]?.[value]) return mappings[key][value];
    if (key === 'threshold') {
      const number = String(value).match(/\d+/)?.[0];
      return number ? `Threshold ${number}` : cellLabel(value);
    }
    return cellLabel(value).replace(/^\d+ /, '');
  }

  function hierarchyRank(key, value) {
    if (key === 'stage') return Number(String(value).match(/^\d+/)?.[0] ?? 999);
    if (key === 'specimen') return ({ALL: 0, CSF: 1, PB: 2})[value] ?? 99;
    if (key === 'threshold') return Number(String(value).match(/\d+/)?.[0] ?? 999);
    if (key === 'strategy') return Number(String(value).match(/^\d+/)?.[0] ?? 999);
    if (key === 'phase') return ({'02_before_preprocess': 0, '04_after_preprocess': 1, '02_data_distribution': 2})[value] ?? 99;
    return 99;
  }

  function sortedGroups(rows) {
    return [...rows].sort((a, b) => {
      for (const key of hierarchyKeys) {
        const difference = hierarchyRank(key, a.dimensions[key]) - hierarchyRank(key, b.dimensions[key]);
        if (difference) return difference;
      }
      return a.id.localeCompare(b.id);
    });
  }

  function hierarchyOptions(key, group) {
    const index = hierarchyKeys.indexOf(key);
    const prefix = hierarchyKeys.slice(0, index);
    const matching = groups.filter(candidate => prefix.every(name =>
      candidate.dimensions[name] === group.dimensions[name]));
    return [...new Set(matching.map(candidate => candidate.dimensions[key]).filter(Boolean))]
      .sort((a, b) => hierarchyRank(key, a) - hierarchyRank(key, b) || a.localeCompare(b));
  }

  function selectHierarchyValue(key, value) {
    const group = currentGroup();
    if (!group) return;
    const index = hierarchyKeys.indexOf(key);
    const constraints = Object.fromEntries(
      hierarchyKeys.slice(0, index).map(name => [name, group.dimensions[name]])
    );
    constraints[key] = value;
    const target = sortedGroups(groups.filter(candidate => Object.entries(constraints)
      .every(([name, expected]) => candidate.dimensions[name] === expected)))[0];
    if (!target) return;
    state.groupId = target.id;
    state.selectedSample = null;
    state.recordId = null;
    window.renderAll();
  }

  function renderHierarchy(group) {
    let visibleIndex = 0;
    const levels = hierarchyKeys.map(key => {
      const options = hierarchyOptions(key, group);
      if (options.length <= 1) return '';
      visibleIndex += 1;
      const [label, description] = hierarchyMeta[key];
      const keyIndex = hierarchyKeys.indexOf(key);
      const prefix = hierarchyKeys.slice(0, keyIndex);
      const buttons = options.map(value => {
        const conditionCount = groups.filter(candidate =>
          prefix.every(name => candidate.dimensions[name] === group.dimensions[name])
          && candidate.dimensions[key] === value).length;
        return `<button data-outlier-dimension="${html(key)}" data-outlier-value="${html(value)}" class="${group.dimensions[key] === value ? 'active' : ''}"><span>${html(hierarchyValueLabel(key, value))}</span><small>${conditionCount} condition${conditionCount === 1 ? '' : 's'}</small></button>`;
      }).join('');
      return `<div class="comparison-level"><div class="comparison-level-label"><span>${visibleIndex}</span><div><b>${html(label)}</b><small>${html(description)}</small></div></div><div class="comparison-switches">${buttons}</div></div>`;
    }).join('');
    return `<section class="comparison-hierarchy outlier-hierarchy" aria-label="Outlier experiment hierarchy">${levels}</section>`;
  }

  function statusCell(cell) {
    if (!cell) return '<td class="outlier-cell na" title="Not measured in this cell-type matrix">NA</td>';
    const flagged = Number(cell.flag) === 1;
    const title = `Distance ${numberText(cell.distance, 3)} · percentile ${percentileText(cell.percentile)}`;
    return `<td class="outlier-cell ${flagged ? 'flagged' : 'normal'}" title="${html(title)}">
      <div class="outlier-cell-status"><b>${flagged ? '1' : '0'}</b><span>${percentileText(cell.percentile)}</span></div>
      <small>d=${numberText(cell.distance, 3)}</small>
    </td>`;
  }

  function renderSampleDetail(sample, group) {
    if (!sample) return '<div class="empty">Select a candidate sample.</div>';
    const review = reviewOf(sample.source_sample_id);
    const cellRows = group.cell_types.map(cellType => {
      const cell = sample.cells[cellType];
      if (!cell) {
        return `<tr><td>${html(cellLabel(cellType))}</td><td><span class="outlier-cell-inline na">NA</span></td><td>NA</td><td>NA</td><td>Not tested</td></tr>`;
      }
      const flagged = Number(cell.flag) === 1;
      const canOpen = cell.record_id && cell.artifact_id;
      return `<tr>
        <td><b>${html(cellLabel(cellType))}</b></td>
        <td><span class="outlier-cell-inline ${flagged ? 'flagged' : 'normal'}">${flagged ? '1' : '0'}</span></td>
        <td title="${html(cell.distance_metric || '')}">${numberText(cell.distance, 3)}</td>
        <td>${percentileText(cell.percentile)}</td>
        <td>${canOpen ? `<button class="outlier-link" data-open-outlier-record="${html(cell.record_id)}" data-open-outlier-artifact="${html(cell.artifact_id)}">Open PCA</button>` : 'Unavailable'}</td>
      </tr>`;
    }).join('');
    return `<article class="outlier-detail outlier-detail-below">
      <div class="outlier-detail-head"><div><h2>${html(sample.source_sample_id)}</h2>
        <div class="sub">${html(sample.cohort || 'cohort NA')} · ${html(sample.diagnosis || 'diagnosis NA')} · ${html(sample.specimen || 'specimen NA')}</div></div>
        <span class="outlier-decision ${html(review.decision)}">${html(decisionLabel(review.decision))}</span></div>
      <div class="outlier-evidence-grid">
        <div><span>Outlier cells</span><b>${numberText(sample.outlier_cell_count, 0)}</b></div>
        <div><span>Weighted count</span><b>${numberText(sample.pca_weighted_count, 2)}</b></div>
        <div><span>Tested cells</span><b>${numberText(sample.tested_cell_count, 0)}</b></div>
        <div><span>Fraction</span><b>${percentileText(100 * sample.outlier_cell_fraction)}</b></div>
      </div>
      <section class="outlier-decision-panel"><div><b>Sample decision</b><div class="decision-row">
        ${[['keep','Keep'],['hold','Hold'],['exclude','Exclude'],['unreviewed','Clear']].map(([value, label]) =>
          `<button data-outlier-decision="${value}" class="${review.decision === value ? 'selected' : ''}">${label}</button>`).join('')}
        </div></div><label><b>Review notes</b><textarea id="outlierNotes" placeholder="Record why this sample should be retained, held, or excluded.">${html(review.notes)}</textarea></label></section>
      <div class="outlier-method-note"><b>Distance는 합산하지 않습니다.</b> Cell type마다 PCA 좌표계와 distance 분포가 다르므로, 비교 가능한 값은 해당 행렬 안의 percentile과 후보 여부(1/0)입니다. Weighted count는 각 후보의 <code>max(0, min(1, (percentile - 95) / 5))</code> 합이며 검토 우선순위용입니다.</div>
      <h3>Cell-type evidence</h3><div class="table-wrap outlier-cell-evidence"><table><thead><tr><th>Cell type</th><th>Outlier</th><th>Distance score</th><th>Percentile</th><th>PCA</th></tr></thead><tbody>${cellRows}</tbody></table></div>
    </article>`;
  }

  function render(container) {
    const group = currentGroup();
    if (!group) {
      container.innerHTML = '<div class="empty">No cell-separated PCA experiment sets were packaged.</div>';
      return;
    }
    const samples = visibleSamples(group);
    if (!state.selectedSample || !samples.some(row => row.source_sample_id === state.selectedSample)) {
      state.selectedSample = samples[0]?.source_sample_id || null;
    }
    const selected = samples.find(row => row.source_sample_id === state.selectedSample) || null;
    const matrixHead = group.cell_types.map(cell => `<th class="outlier-cell-head">${html(cellLabel(cell))}</th>`).join('');
    const matrixRows = samples.map((sample, index) => {
      const review = reviewOf(sample.source_sample_id);
      return `<tr data-outlier-sample="${html(sample.source_sample_id)}" class="${sample.source_sample_id === state.selectedSample ? 'selected' : ''}">
        <td>${index + 1}</td><td><b>${html(sample.source_sample_id)}</b><div class="table-secondary">${html(sample.diagnosis || 'diagnosis NA')} · ${html(sample.cohort || 'cohort NA')}</div></td>
        ${group.cell_types.map(cell => statusCell(sample.cells[cell])).join('')}
        <td>${numberText(sample.outlier_cell_count, 0)}</td><td>${numberText(sample.pca_weighted_count, 2)}</td><td>${numberText(sample.tested_cell_count, 0)}</td><td>${percentileText(100 * sample.outlier_cell_fraction)}</td>
        <td><span class="outlier-decision ${html(review.decision)}">${html(decisionLabel(review.decision))}</span></td>
      </tr>`;
    }).join('');
    const dimensions = group.dimensions || {};
    const context = [dimensions.stage, dimensions.threshold, dimensions.strategy, dimensions.specimen, dimensions.phase]
      .filter(Boolean).join(' · ');
    const recordBanner = state.recordId
      ? `<div class="outlier-context">Opened from an individual PCA run. The matched experiment set is selected.<button id="outlierClearRecord">Keep selection and clear shortcut</button></div>` : '';
    container.innerHTML = `<div class="comparison-head"><div><h2>Outlier candidates by cell type</h2>
      <div class="sub">One table per matched experiment condition; saved PCA candidate flags are shown without recalculation.</div></div>
      <div class="comparison-counts">${numberText(samples.length, 0)} candidates · ${numberText(group.cell_types.length, 0)} cell types</div></div>
      ${recordBanner}
      ${renderHierarchy(group)}
      <section class="outlier-controls outlier-filter-controls">
        <label><span>Sample search</span><input id="outlierSearch" type="search" placeholder="SampleID, cohort, diagnosis" value="${html(state.query)}"></label>
        <label><span>Decision</span><select id="outlierDecision"><option value="">All decisions</option>${['unreviewed','keep','hold','exclude'].map(value => `<option value="${value}" ${state.decision === value ? 'selected' : ''}>${decisionLabel(value)}</option>`).join('')}</select></label>
        <label><span>Ordering</span><select id="outlierSort"><option value="weighted" ${state.sort === 'weighted' ? 'selected' : ''}>Highest weighted count</option><option value="count" ${state.sort === 'count' ? 'selected' : ''}>Most outlier cells</option><option value="fraction" ${state.sort === 'fraction' ? 'selected' : ''}>Highest fraction</option><option value="sample" ${state.sort === 'sample' ? 'selected' : ''}>Sample ID</option></select></label>
      </section>
      <div class="outlier-downloads"><a href="${html(data.downloads.samples)}" download>Sample counts CSV</a><a href="${html(data.downloads.occurrences)}" download>Occurrence evidence CSV</a><a href="${html(data.downloads.experiments)}" download>Experiment counts CSV</a><button id="outlierExportReview">Export review decisions</button></div>
      <p class="outlier-matrix-note"><b>${html(context)}</b> · 각 cell은 후보 여부(1/0), 같은 PCA 안의 percentile, raw distance(d)를 함께 표시합니다. NA는 해당 cell-type PCA에 샘플이 없다는 뜻입니다. Cell type 간 심각도 비교에는 percentile을 사용하고 raw distance는 해당 PCA 내부에서만 해석합니다.</p>
      <div class="table-wrap outlier-candidate-table outlier-matrix-table"><table><thead><tr><th>#</th><th>SampleID</th>${matrixHead}<th>Outlier cells</th><th>Weighted count</th><th>Tested cells</th><th>Fraction</th><th>Decision</th></tr></thead><tbody>${matrixRows || '<tr><td colspan="99">No candidate matches the current filters.</td></tr>'}</tbody></table></div>
      ${renderSampleDetail(selected, group)}`;

    document.querySelectorAll('[data-outlier-dimension]').forEach(node => node.onclick = () =>
      selectHierarchyValue(node.dataset.outlierDimension, node.dataset.outlierValue));
    document.getElementById('outlierDecision').onchange = event => { state.decision = event.target.value; state.selectedSample = null; window.renderAll(); };
    document.getElementById('outlierSort').onchange = event => { state.sort = event.target.value; window.renderAll(); };
    document.getElementById('outlierSearch').oninput = event => {
      state.query = event.target.value; clearTimeout(searchTimer);
      searchTimer = setTimeout(() => {
        render(container); renderKpis();
        const input = document.getElementById('outlierSearch'); input?.focus();
        input?.setSelectionRange(input.value.length, input.value.length);
      }, 160);
    };
    document.querySelectorAll('[data-outlier-sample]').forEach(node => node.onclick = () => {
      state.selectedSample = node.dataset.outlierSample; render(container);
    });
    document.querySelectorAll('[data-outlier-decision]').forEach(node => node.onclick = () => {
      if (selected) persistReview(selected.source_sample_id, {decision: node.dataset.outlierDecision});
    });
    const notes = document.getElementById('outlierNotes');
    if (notes && selected) notes.onchange = event => persistReview(selected.source_sample_id, {notes: event.target.value});
    document.querySelectorAll('[data-open-outlier-record]').forEach(node => node.onclick = () =>
      window.openOutlierRecord?.(node.dataset.openOutlierRecord, node.dataset.openOutlierArtifact));
    document.getElementById('outlierExportReview').onclick = exportReviews;
    const clear = document.getElementById('outlierClearRecord');
    if (clear) clear.onclick = () => { state.recordId = null; render(container); };
  }

  function csvCell(value) {
    const text = String(value ?? '');
    return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  }

  function exportReviews() {
    const header = ['source_sample_id', 'decision', 'notes', 'updated_at'];
    const rows = Object.entries(reviews)
      .filter(([, value]) => value.decision !== 'unreviewed' || value.notes)
      .map(([sampleId, value]) => [sampleId, value.decision, value.notes, value.updated_at]);
    const blob = new Blob([[header, ...rows].map(row => row.map(csvCell).join(',')).join('\n')], {type: 'text/csv;charset=utf-8'});
    const anchor = document.createElement('a');
    anchor.href = URL.createObjectURL(blob); anchor.download = `outlier_review_${new Date().toISOString().slice(0, 10)}.csv`;
    anchor.click(); URL.revokeObjectURL(anchor.href);
  }

  function groupForRecord(recordId) {
    return groups.find(group => group.candidate_samples.some(sample =>
      Object.values(sample.cells).some(cell => cell.record_id === recordId)));
  }

  function openForRecord(recordId) {
    const group = groupForRecord(recordId);
    if (group) state.groupId = group.id;
    state.recordId = recordId; state.selectedSample = null;
    document.getElementById('viewMode').value = 'outliers'; window.renderAll();
  }

  function attachRunShortcut(recordId) {
    const actions = document.querySelector('#detail .detail-head .actions');
    if (!actions || !recordId) return;
    const group = groupForRecord(recordId);
    if (!group) return;
    const count = group.candidate_samples.filter(sample =>
      Object.values(sample.cells).some(cell => cell.record_id === recordId && Number(cell.flag) === 1)).length;
    const button = document.createElement('button');
    button.type = 'button'; button.textContent = `Outlier candidates (${count})`;
    button.onclick = () => openForRecord(recordId); actions.prepend(button);
  }

  window.OutlierReview = {render, renderKpis, openForRecord, attachRunShortcut, data};
})();
