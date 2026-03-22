const CARD_METRICS = [
  { key: 'num_chains', label: 'Chains', format: 'int' },
  { key: 'num_residues', label: 'Residues', format: 'int' },
  { key: 'num_sequence_unique_chains', label: 'Unique chains', format: 'int' },
  { key: 'frac_alpha', label: 'Alpha', format: 'percent' },
  { key: 'frac_beta', label: 'Beta', format: 'percent' },
  { key: 'frac_coil', label: 'Coil', format: 'percent' },
  { key: 'num_pores', label: 'Pores', format: 'int' },
  { key: 'largest_pore_volume_a3', label: 'Pore vol.', format: 'float0' },
  { key: 'largest_pore_length_a', label: 'Pore len.', format: 'float1' },
  { key: 'largest_pore_min_diameter_a', label: 'Pore d-min', format: 'float1' },
  { key: 'largest_pore_max_diameter_a', label: 'Pore d-max', format: 'float1' },
  { key: 'largest_pore_circularity', label: 'Pore circ.', format: 'float2' },
  { key: 'largest_pore_uniformity', label: 'Pore unif.', format: 'float2' },
  { key: 'num_pockets', label: 'Pockets', format: 'int' },
  { key: 'largest_pocket_volume_a3', label: 'Pocket vol.', format: 'float0' },
  { key: 'largest_pocket_length_a', label: 'Pocket len.', format: 'float1' },
  { key: 'largest_pocket_min_diameter_a', label: 'Pocket d-min', format: 'float1' },
  { key: 'largest_pocket_max_diameter_a', label: 'Pocket d-max', format: 'float1' },
  { key: 'largest_pocket_circularity', label: 'Pocket circ.', format: 'float2' },
  { key: 'largest_pocket_uniformity', label: 'Pocket unif.', format: 'float2' },
  { key: 'num_cavities', label: 'Cavities', format: 'int' },
  { key: 'largest_cavity_volume_a3', label: 'Cavity vol.', format: 'float0' },
  { key: 'largest_cavity_length_a', label: 'Cavity len.', format: 'float1' },
  { key: 'largest_cavity_min_diameter_a', label: 'Cavity d-min', format: 'float1' },
  { key: 'largest_cavity_max_diameter_a', label: 'Cavity d-max', format: 'float1' },
  { key: 'largest_cavity_circularity', label: 'Cavity circ.', format: 'float2' },
  { key: 'largest_cavity_uniformity', label: 'Cavity unif.', format: 'float2' },
  { key: 'num_hubs', label: 'Hubs', format: 'int' },
  { key: 'largest_hub_volume_a3', label: 'Hub vol.', format: 'float0' },
  { key: 'largest_hub_length_a', label: 'Hub len.', format: 'float1' },
  { key: 'largest_hub_min_diameter_a', label: 'Hub d-min', format: 'float1' },
  { key: 'largest_hub_max_diameter_a', label: 'Hub d-max', format: 'float1' },
  { key: 'largest_hub_circularity', label: 'Hub circ.', format: 'float2' },
  { key: 'largest_hub_uniformity', label: 'Hub unif.', format: 'float2' },
];

const FILTER_PRESETS_KEY = 'volumizer_filter_presets';
const DISPLAY_PRESETS_KEY = 'volumizer_display_presets';

const state = {
  currentOffset: 0,
  totalCount: 0,
  currentLimit: 24,
  currentViewer: null,
};

const elements = {
  filtersForm: document.getElementById('filters-form'),
  runSelect: document.getElementById('run-id'),
  limitSelect: document.getElementById('limit'),
  sortBySelect: document.getElementById('sort-by'),
  sortDirSelect: document.getElementById('sort-dir'),
  prevPageButton: document.getElementById('prev-page'),
  nextPageButton: document.getElementById('next-page'),
  resetButton: document.getElementById('reset-filters'),
  loadingState: document.getElementById('loading-state'),
  emptyState: document.getElementById('empty-state'),
  resultsGrid: document.getElementById('results-grid'),
  dbPill: document.getElementById('db-pill'),
  resultPill: document.getElementById('result-pill'),
  pageStatus: document.getElementById('page-status'),
  detailDialog: document.getElementById('detail-dialog'),
  closeDetailButton: document.getElementById('close-detail'),
  detailTitle: document.getElementById('detail-title'),
  detailLinks: document.getElementById('detail-links'),
  detailMetrics: document.getElementById('detail-metrics'),
  detailVolumes: document.getElementById('detail-volumes'),
  viewerHost: document.getElementById('viewer-host'),
  displayPropsToggle: document.getElementById('display-props-toggle'),
  displayPropsPopup: document.getElementById('display-props-popup'),
  filterPresetSelect: document.getElementById('filter-preset-select'),
  filterPresetLoad: document.getElementById('filter-preset-load'),
  filterPresetSave: document.getElementById('filter-preset-save'),
  filterPresetDelete: document.getElementById('filter-preset-delete'),
  displayPresetSelect: document.getElementById('display-preset-select'),
  displayPresetLoad: document.getElementById('display-preset-load'),
  displayPresetSave: document.getElementById('display-preset-save'),
  displayPresetDelete: document.getElementById('display-preset-delete'),
};

function numericInputValue(name) {
  const field = elements.filtersForm.elements.namedItem(name);
  const value = field ? String(field.value).trim() : '';
  if (value === '') return null;
  return value;
}

function buildSearchParams() {
  const params = new URLSearchParams();
  const runId = String(elements.runSelect.value || '').trim();
  if (runId) params.set('run_id', runId);

  const numericNames = [
    'pore_volume_min',
    'pore_length_min',
    'pore_dmin_min',
    'pore_dmax_min',
    'pocket_volume_min',
    'pocket_length_min',
    'pocket_dmin_min',
    'pocket_dmax_min',
    'cavity_volume_min',
    'cavity_length_min',
    'cavity_dmin_min',
    'cavity_dmax_min',
    'hub_volume_min',
    'hub_length_min',
    'hub_dmin_min',
    'hub_dmax_min',
    'chains_min',
    'residues_min',
    'seq_unique_chains_min',
    'frac_alpha_min',
    'frac_alpha_max',
    'frac_beta_min',
    'frac_beta_max',
    'frac_coil_min',
    'frac_coil_max',
    'pore_circularity_min',
    'pore_uniformity_min',
    'pocket_circularity_min',
    'pocket_uniformity_min',
    'cavity_circularity_min',
    'cavity_uniformity_min',
    'hub_circularity_min',
    'hub_uniformity_min',
  ];
  for (const name of numericNames) {
    const value = numericInputValue(name);
    if (value !== null) params.set(name, value);
  }

  params.set('sort_by', elements.sortBySelect.value);
  params.set('sort_dir', elements.sortDirSelect.value);
  params.set('limit', String(state.currentLimit));
  params.set('offset', String(state.currentOffset));
  return params;
}

async function fetchJson(url) {
  const response = await fetch(url, { headers: { Accept: 'application/json' } });
  if (!response.ok) {
    let detail = `${response.status} ${response.statusText}`;
    try {
      const payload = await response.json();
      if (payload && payload.detail) detail = String(payload.detail);
    } catch {
      // Ignore JSON parse failures for plain-text responses.
    }
    throw new Error(detail);
  }
  return response.json();
}

async function loadHealth() {
  const health = await fetchJson('/api/health');
  elements.dbPill.textContent = health.db_exists ? health.db.split('/').slice(-2).join('/') : 'Missing DB';
}

async function loadRuns() {
  const payload = await fetchJson('/api/runs');
  const previousValue = elements.runSelect.value;
  elements.runSelect.innerHTML = '<option value="">All runs</option>';

  for (const run of payload.runs) {
    const option = document.createElement('option');
    option.value = run.run_id;
    option.textContent = `${run.run_id} (${run.structure_count})`;
    elements.runSelect.append(option);
  }

  if ([...elements.runSelect.options].some((option) => option.value === previousValue)) {
    elements.runSelect.value = previousValue;
  }
}

function renderThumbnail(url, axis) {
  const frame = document.createElement('div');
  frame.className = 'thumb-frame';

  if (url) {
    const image = document.createElement('img');
    image.src = url;
    image.alt = `${axis.toUpperCase()} thumbnail`;
    frame.append(image);
  } else {
    const placeholder = document.createElement('div');
    placeholder.className = 'thumb-placeholder';
    placeholder.textContent = axis.toUpperCase();
    frame.append(placeholder);
  }

  return frame;
}

function renderMetric(label, value) {
  const item = document.createElement('div');
  item.className = 'metric-pill';
  const labelNode = document.createElement('span');
  labelNode.className = 'metric-label';
  labelNode.textContent = label;
  const valueNode = document.createElement('span');
  valueNode.className = 'metric-value';
  valueNode.textContent = value;
  item.append(labelNode, valueNode);
  return item;
}

function prettyNumber(value, digits = 1) {
  if (value === null || value === undefined || value === '') return '—';
  const numericValue = Number(value);
  if (!Number.isFinite(numericValue)) return String(value);
  return numericValue.toLocaleString(undefined, {
    maximumFractionDigits: digits,
  });
}

function prettyPercent(value) {
  if (value === null || value === undefined || value === '') return '—';
  return `${(Number(value) * 100).toFixed(0)}%`;
}

function getVisibleMetricKeys() {
  const keys = new Set();
  for (const checkbox of elements.displayPropsPopup.querySelectorAll('input[data-metric]')) {
    if (checkbox.checked) keys.add(checkbox.dataset.metric);
  }
  return keys;
}

function formatMetricValue(value, format) {
  if (format === 'int') return prettyNumber(value, 0);
  if (format === 'percent') return prettyPercent(value);
  if (format === 'float0') return prettyNumber(value, 0);
  if (format === 'float1') return prettyNumber(value, 1);
  if (format === 'float2') return prettyNumber(value, 2);
  return prettyNumber(value, 1);
}

function renderResults(rows) {
  elements.resultsGrid.innerHTML = '';
  for (const row of rows) {
    const card = document.createElement('article');
    card.className = 'result-card';

    const thumbStrip = document.createElement('div');
    thumbStrip.className = 'thumb-strip';
    thumbStrip.append(
      renderThumbnail(row.urls.thumbnail_x, 'x'),
      renderThumbnail(row.urls.thumbnail_y, 'y'),
      renderThumbnail(row.urls.thumbnail_z, 'z'),
    );

    const body = document.createElement('div');
    body.className = 'card-body';

    const titleRow = document.createElement('div');
    titleRow.className = 'card-title-row';
    const title = document.createElement('h3');
    title.textContent = row.source_label;
    const runTag = document.createElement('span');
    runTag.className = 'run-tag';
    runTag.textContent = row.run_id;
    titleRow.append(title, runTag);

    const subtitle = document.createElement('p');
    subtitle.className = 'card-subtitle';
    subtitle.textContent = row.pdb_id ? `PDB ${row.pdb_id}` : 'Local structure';

    const metrics = document.createElement('div');
    metrics.className = 'metric-grid';
    const visibleKeys = getVisibleMetricKeys();
    for (const metric of CARD_METRICS) {
      if (visibleKeys.has(metric.key)) {
        metrics.append(renderMetric(metric.label, formatMetricValue(row[metric.key], metric.format)));
      }
    }

    const actions = document.createElement('div');
    actions.className = 'card-actions';
    const openButton = document.createElement('button');
    openButton.type = 'button';
    openButton.className = 'primary-button';
    openButton.textContent = 'Open Detail';
    openButton.addEventListener('click', () => openDetail(row.structure_id));
    actions.append(openButton);

    body.append(titleRow, subtitle, metrics, actions);
    card.append(thumbStrip, body);
    elements.resultsGrid.append(card);
  }
}

function updatePaging() {
  const pageNumber = Math.floor(state.currentOffset / state.currentLimit) + 1;
  const pageCount = Math.max(1, Math.ceil(state.totalCount / state.currentLimit));
  elements.pageStatus.textContent = `Page ${pageNumber} of ${pageCount}`;
  elements.prevPageButton.disabled = state.currentOffset <= 0;
  elements.nextPageButton.disabled = state.currentOffset + state.currentLimit >= state.totalCount;
}

async function search() {
  elements.loadingState.hidden = false;
  elements.emptyState.hidden = true;

  try {
    const params = buildSearchParams();
    const payload = await fetchJson(`/api/hits?${params.toString()}`);
    state.totalCount = payload.total_count;
    elements.resultPill.textContent = String(payload.total_count);
    renderResults(payload.rows);
    elements.emptyState.hidden = payload.rows.length > 0;
    updatePaging();
  } catch (error) {
    elements.resultsGrid.innerHTML = '';
    elements.emptyState.hidden = false;
    elements.emptyState.textContent = `Unable to load hits: ${error.message}`;
    state.totalCount = 0;
    elements.resultPill.textContent = '0';
    updatePaging();
  } finally {
    elements.loadingState.hidden = true;
  }
}

function buildDetailMetric(label, value) {
  const wrapper = document.createElement('div');
  wrapper.className = 'detail-metric';
  const labelNode = document.createElement('span');
  labelNode.className = 'detail-metric-label';
  labelNode.textContent = label;
  const valueNode = document.createElement('strong');
  valueNode.textContent = value;
  wrapper.append(labelNode, valueNode);
  return wrapper;
}

function renderDetail(detail, viewerData) {
  elements.detailTitle.textContent = detail.source_label;
  elements.detailLinks.innerHTML = '';

  const linkItems = [
    ['Annotated CIF', viewerData.structure_url],
    ['Annotation JSON', viewerData.annotation_url],
  ];
  for (const [label, href] of linkItems) {
    if (!href) continue;
    const link = document.createElement('a');
    link.href = href;
    link.target = '_blank';
    link.rel = 'noreferrer';
    link.textContent = label;
    elements.detailLinks.append(link);
  }

  elements.detailMetrics.innerHTML = '';
  elements.detailMetrics.append(
    buildDetailMetric('Run', detail.run_id),
    buildDetailMetric('PDB', detail.pdb_id || '—'),
    buildDetailMetric('Pores', prettyNumber(detail.num_pores, 0)),
    buildDetailMetric('Pockets', prettyNumber(detail.num_pockets, 0)),
    buildDetailMetric('Cavities', prettyNumber(detail.num_cavities, 0)),
    buildDetailMetric('Hubs', prettyNumber(detail.num_hubs, 0)),
    buildDetailMetric('Chains', prettyNumber(detail.num_chains, 0)),
    buildDetailMetric('Residues', prettyNumber(detail.num_residues, 0)),
    buildDetailMetric('Unique chains', prettyNumber(detail.num_sequence_unique_chains, 0)),
    buildDetailMetric('Alpha', prettyPercent(detail.frac_alpha)),
    buildDetailMetric('Beta', prettyPercent(detail.frac_beta)),
    buildDetailMetric('Coil', prettyPercent(detail.frac_coil)),
    buildDetailMetric('Pore circ.', prettyNumber(detail.largest_pore_circularity, 2)),
    buildDetailMetric('Pore unif.', prettyNumber(detail.largest_pore_uniformity, 2)),
    buildDetailMetric('Pocket circ.', prettyNumber(detail.largest_pocket_circularity, 2)),
    buildDetailMetric('Pocket unif.', prettyNumber(detail.largest_pocket_uniformity, 2)),
  );

  elements.detailVolumes.innerHTML = '';
  for (const volume of detail.volumes) {
    const row = document.createElement('tr');
    row.innerHTML = [
      `<td>${volume.kind}</td>`,
      `<td>${Number(volume.rank_in_kind) + 1}</td>`,
      `<td>${prettyNumber(volume.volume_a3, 0)}</td>`,
      `<td>${prettyNumber(volume.length_a)}</td>`,
      `<td>${prettyNumber(volume.min_diameter_a)}</td>`,
      `<td>${prettyNumber(volume.max_diameter_a)}</td>`,
      `<td>${prettyNumber(volume.cross_section_circularity, 2)}</td>`,
      `<td>${prettyNumber(volume.cross_section_uniformity, 2)}</td>`,
    ].join('');
    elements.detailVolumes.append(row);
  }
}

async function mountViewer(viewerData) {
  if (!viewerData.structure_url || !viewerData.structure_format) {
    elements.viewerHost.innerHTML = '<div class="viewer-placeholder">Annotated structure is unavailable for this hit.</div>';
    return;
  }

  if (!window.molstar || !window.molstar.Viewer) {
    elements.viewerHost.innerHTML = '<div class="viewer-placeholder">Mol* failed to load. Open the annotated CIF directly instead.</div>';
    return;
  }

  if (state.currentViewer && state.currentViewer.plugin && typeof state.currentViewer.plugin.dispose === 'function') {
    state.currentViewer.plugin.dispose();
    state.currentViewer = null;
  }

  elements.viewerHost.innerHTML = '<div id="molstar-app" class="molstar-app"></div>';

  const viewer = await window.molstar.Viewer.create('molstar-app', {
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLog: false,
    layoutShowLeftPanel: false,
    viewportShowExpand: false,
    viewportShowSelectionMode: false,
    viewportShowAnimation: false,
    collapseLeftPanel: true,
    pluginStateServer: false,
  });

  state.currentViewer = viewer;
  await viewer.loadStructureFromUrl(
    viewerData.structure_url,
    viewerData.structure_format,
    viewerData.structure_format === 'bcif',
  );
}

async function openDetail(structureId) {
  try {
    const [detail, viewerData] = await Promise.all([
      fetchJson(`/api/hits/${structureId}`),
      fetchJson(`/api/hits/${structureId}/viewer-data`),
    ]);
    renderDetail(detail, viewerData);
    if (!elements.detailDialog.open) {
      elements.detailDialog.showModal();
    }
    await mountViewer(viewerData);
  } catch (error) {
    elements.viewerHost.innerHTML = `<div class="viewer-placeholder">Unable to load detail: ${error.message}</div>`;
    if (!elements.detailDialog.open) {
      elements.detailDialog.showModal();
    }
  }
}

function closeDetail() {
  if (elements.detailDialog.open) {
    elements.detailDialog.close();
  }
}

/* ---- Presets ---- */

function loadPresetsMap(storageKey) {
  try {
    return JSON.parse(localStorage.getItem(storageKey)) || {};
  } catch {
    return {};
  }
}

function savePresetsMap(storageKey, map) {
  localStorage.setItem(storageKey, JSON.stringify(map));
}

function populatePresetSelect(selectEl, storageKey) {
  const map = loadPresetsMap(storageKey);
  const previousValue = selectEl.value;
  selectEl.innerHTML = '<option value="">Presets\u2026</option>';
  for (const name of Object.keys(map).sort()) {
    const opt = document.createElement('option');
    opt.value = name;
    opt.textContent = name;
    selectEl.append(opt);
  }
  if ([...selectEl.options].some((o) => o.value === previousValue)) {
    selectEl.value = previousValue;
  }
}

function captureFilterState() {
  const data = {};
  for (const el of elements.filtersForm.elements) {
    if (el.name) data[el.name] = el.value;
  }
  return data;
}

function applyFilterState(data) {
  for (const el of elements.filtersForm.elements) {
    if (el.name && data[el.name] !== undefined) {
      el.value = data[el.name];
    }
  }
  state.currentLimit = Number(elements.limitSelect.value);
  state.currentOffset = 0;
}

function captureDisplayState() {
  const keys = [];
  for (const cb of elements.displayPropsPopup.querySelectorAll('input[data-metric]')) {
    if (cb.checked) keys.push(cb.dataset.metric);
  }
  return keys;
}

function applyDisplayState(keys) {
  const set = new Set(keys);
  for (const cb of elements.displayPropsPopup.querySelectorAll('input[data-metric]')) {
    cb.checked = set.has(cb.dataset.metric);
  }
}

function promptPresetName(action) {
  const name = prompt(`Enter a name to ${action}:`);
  if (!name || !name.trim()) return null;
  return name.trim();
}

function wirePresets() {
  // Filter presets
  populatePresetSelect(elements.filterPresetSelect, FILTER_PRESETS_KEY);

  elements.filterPresetSave.addEventListener('click', () => {
    const name = promptPresetName('save this filter preset');
    if (!name) return;
    const map = loadPresetsMap(FILTER_PRESETS_KEY);
    if (map[name] && !confirm(`Overwrite existing preset "${name}"?`)) return;
    map[name] = captureFilterState();
    savePresetsMap(FILTER_PRESETS_KEY, map);
    populatePresetSelect(elements.filterPresetSelect, FILTER_PRESETS_KEY);
    elements.filterPresetSelect.value = name;
  });

  elements.filterPresetLoad.addEventListener('click', () => {
    const name = elements.filterPresetSelect.value;
    if (!name) return;
    const map = loadPresetsMap(FILTER_PRESETS_KEY);
    if (!map[name]) return;
    applyFilterState(map[name]);
    search();
  });

  elements.filterPresetDelete.addEventListener('click', () => {
    const name = elements.filterPresetSelect.value;
    if (!name) return;
    if (!confirm(`Delete preset "${name}"?`)) return;
    const map = loadPresetsMap(FILTER_PRESETS_KEY);
    delete map[name];
    savePresetsMap(FILTER_PRESETS_KEY, map);
    populatePresetSelect(elements.filterPresetSelect, FILTER_PRESETS_KEY);
  });

  // Display property presets
  populatePresetSelect(elements.displayPresetSelect, DISPLAY_PRESETS_KEY);

  elements.displayPresetSave.addEventListener('click', () => {
    const name = promptPresetName('save this display preset');
    if (!name) return;
    const map = loadPresetsMap(DISPLAY_PRESETS_KEY);
    if (map[name] && !confirm(`Overwrite existing preset "${name}"?`)) return;
    map[name] = captureDisplayState();
    savePresetsMap(DISPLAY_PRESETS_KEY, map);
    populatePresetSelect(elements.displayPresetSelect, DISPLAY_PRESETS_KEY);
    elements.displayPresetSelect.value = name;
  });

  elements.displayPresetLoad.addEventListener('click', () => {
    const name = elements.displayPresetSelect.value;
    if (!name) return;
    const map = loadPresetsMap(DISPLAY_PRESETS_KEY);
    if (!map[name]) return;
    applyDisplayState(map[name]);
    search();
  });

  elements.displayPresetDelete.addEventListener('click', () => {
    const name = elements.displayPresetSelect.value;
    if (!name) return;
    if (!confirm(`Delete preset "${name}"?`)) return;
    const map = loadPresetsMap(DISPLAY_PRESETS_KEY);
    delete map[name];
    savePresetsMap(DISPLAY_PRESETS_KEY, map);
    populatePresetSelect(elements.displayPresetSelect, DISPLAY_PRESETS_KEY);
  });
}

function wireEvents() {
  elements.filtersForm.addEventListener('submit', async (event) => {
    event.preventDefault();
    state.currentOffset = 0;
    state.currentLimit = Number(elements.limitSelect.value);
    await search();
  });

  elements.resetButton.addEventListener('click', async () => {
    elements.sortBySelect.value = 'largest_pore_volume';
    elements.sortDirSelect.value = 'desc';
    elements.limitSelect.value = '24';
    state.currentLimit = 24;
    state.currentOffset = 0;
    window.setTimeout(search, 0);
  });

  elements.prevPageButton.addEventListener('click', async () => {
    state.currentOffset = Math.max(0, state.currentOffset - state.currentLimit);
    await search();
  });

  elements.nextPageButton.addEventListener('click', async () => {
    state.currentOffset += state.currentLimit;
    await search();
  });

  elements.displayPropsToggle.addEventListener('click', () => {
    elements.displayPropsPopup.hidden = !elements.displayPropsPopup.hidden;
  });

  for (const checkbox of elements.displayPropsPopup.querySelectorAll('input[data-metric]')) {
    checkbox.addEventListener('change', () => search());
  }

  elements.closeDetailButton.addEventListener('click', closeDetail);
  elements.detailDialog.addEventListener('click', (event) => {
    const rect = elements.detailDialog.getBoundingClientRect();
    const insideDialog =
      rect.top <= event.clientY &&
      event.clientY <= rect.top + rect.height &&
      rect.left <= event.clientX &&
      event.clientX <= rect.left + rect.width;
    if (!insideDialog) {
      closeDetail();
    }
  });
}

async function init() {
  wireEvents();
  wirePresets();
  await loadHealth();
  await loadRuns();
  await search();
}

window.addEventListener('DOMContentLoaded', () => {
  init().catch((error) => {
    elements.loadingState.hidden = true;
    elements.emptyState.hidden = false;
    elements.emptyState.textContent = `Startup failed: ${error.message}`;
  });
});
