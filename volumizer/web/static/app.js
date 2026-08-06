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
];

const FILTER_PRESETS_KEY = 'volumizer_filter_presets';
const DISPLAY_PRESETS_KEY = 'volumizer_display_presets';
const ACTIVE_FILTER_PRESET_KEY = 'volumizer_active_filter_preset';
const FILTER_PRESET_NONE_LABEL = 'None';
const BUILTIN_FILTER_PRESET_PREFIX = '__builtin__:';
const DISPLAY_PRESET_DEFAULT_LABEL = 'Default';
const BUILTIN_DISPLAY_PRESET_PREFIX = '__builtin-display__:';
const NON_FILTER_FORM_FIELD_NAMES = new Set(['sort_by', 'sort_dir', 'limit']);
const DEFAULT_UPLOAD_RESOLUTION = '3.0';
const DEFAULT_ASSEMBLY_POLICIES = Object.freeze(['biological', 'asymmetric', 'auto']);
const UPLOAD_POLL_INTERVAL_MS = 1000;
const SORT_BY_ALIASES = Object.freeze({
  largest_pore_volume_a3: 'largest_pore_volume',
  largest_pore_length_a: 'largest_pore_length',
  largest_pore_min_diameter_a: 'largest_pore_min_diameter',
  largest_pore_max_diameter_a: 'largest_pore_max_diameter',
  largest_pocket_volume_a3: 'largest_pocket_volume',
  largest_pocket_length_a: 'largest_pocket_length',
  largest_pocket_min_diameter_a: 'largest_pocket_min_diameter',
  largest_pocket_max_diameter_a: 'largest_pocket_max_diameter',
  largest_cavity_volume_a3: 'largest_cavity_volume',
  largest_cavity_length_a: 'largest_cavity_length',
  largest_cavity_min_diameter_a: 'largest_cavity_min_diameter',
  largest_cavity_max_diameter_a: 'largest_cavity_max_diameter',
});
const DEFAULT_DISPLAY_METRICS = Object.freeze([
  'num_chains',
  'num_residues',
  'num_sequence_unique_chains',
]);
const BUILTIN_FILTER_PRESETS = Object.freeze({
  Pores: Object.freeze({
    seq_unique_chains_max: '3',
    frac_coil_max: '0.6',
    pore_volume_min: '10000',
    pore_dmin_min: '35',
    pore_circularity_min: '0.8',
    pore_uniformity_min: '0.5',
  }),
  Pockets: Object.freeze({
    seq_unique_chains_max: '3',
    frac_coil_max: '0.6',
  }),
  Cavities: Object.freeze({
    seq_unique_chains_max: '3',
    frac_coil_max: '0.6',
  }),
});
const BUILTIN_DISPLAY_PRESETS = Object.freeze({
  Pores: Object.freeze([
    ...DEFAULT_DISPLAY_METRICS,
    'largest_pore_volume_a3',
    'largest_pore_min_diameter_a',
    'largest_pore_max_diameter_a',
    'largest_pore_uniformity',
    'largest_pore_circularity',
  ]),
  Pockets: Object.freeze([
    ...DEFAULT_DISPLAY_METRICS,
    'largest_pocket_volume_a3',
    'largest_pocket_max_diameter_a',
    'largest_pocket_length_a',
    'largest_pocket_circularity',
  ]),
  Cavities: Object.freeze([
    ...DEFAULT_DISPLAY_METRICS,
    'largest_cavity_volume_a3',
    'largest_cavity_max_diameter_a',
    'largest_cavity_length_a',
    'largest_cavity_circularity',
    'largest_cavity_uniformity',
  ]),
});

const NON_PROTEIN_COMP_IDS = ['HUB', 'POR', 'POK', 'CAV', 'OCC'];
const VOLUME_STYLES = {
  POR: {
    label: 'Pores',
    color: 0xDD8833,      // orange
    mouthColor: 0xA85618, // dark orange
    mouthLabel: 'Pore Mouths',
  },
  POK: {
    label: 'Pockets',
    color: 0x3366CC,      // blue
    mouthColor: 0x1F3F80, // dark blue
    mouthLabel: 'Pocket Mouths',
  },
  CAV: {
    label: 'Cavities',
    color: 0xCC33CC,      // magenta
    mouthColor: 0x8A238A, // dark magenta
    mouthLabel: 'Cavity Mouths',
  },
};

const PROTEIN_PALETTE = [
  0x2E8B57, 0x808080, 0x3CB371, 0xA0A0A0,
  0x228B22, 0x6B6B6B, 0x006400, 0xB5B5B5,
  0x4CAF50, 0x959595, 0x66BB6A, 0x7A7A7A,
];

// MolScript text expression selecting displayed volumes and optional shell flag.
function volumeScript(compId, variant = 'all') {
  let atomTest = '';
  if (variant === 'core') {
    atomTest = ' :atom-test (<= atom.B_iso_or_equiv 0)';
  } else if (variant === 'mouth') {
    atomTest = ' :atom-test (> atom.B_iso_or_equiv 0)';
  }

  const style = VOLUME_STYLES[compId];
  const label = variant === 'mouth'
    ? style?.mouthLabel || `${compId} Mouths`
    : style?.label || compId;

  return {
    type: {
      name: 'script',
      params: {
        language: 'mol-script',
        expression: `(sel.atom.atom-groups :residue-test (= atom.label_comp_id ${compId})${atomTest})`,
      },
    },
    nullIfEmpty: true,
    label,
  };
}

const state = {
  currentView: 'landing',
  currentOffset: 0,
  totalCount: 0,
  currentLimit: 24,
  currentViewer: null,
  health: null,
  uploadPollId: null,
  uploadRoutedJobId: null,
  suppressDetailCloseRoute: false,
};

const elements = {
  landingView: document.getElementById('landing-view'),
  galleryView: document.getElementById('gallery-view'),
  homeViewButton: document.getElementById('home-view-button'),
  galleryViewButton: document.getElementById('gallery-view-button'),
  uploadFocusButton: document.getElementById('upload-focus-button'),
  browseGalleryButton: document.getElementById('browse-gallery-button'),
  uploadForm: document.getElementById('upload-form'),
  uploadFileInput: document.getElementById('upload-file'),
  uploadResolutionInput: document.getElementById('upload-resolution'),
  uploadAssemblyPolicySelect: document.getElementById('upload-assembly-policy'),
  uploadSubmitButton: document.getElementById('upload-submit-button'),
  uploadResetButton: document.getElementById('upload-reset-button'),
  uploadDisabledState: document.getElementById('upload-disabled-state'),
  uploadStatusPanel: document.getElementById('upload-status-panel'),
  uploadStatusTitle: document.getElementById('upload-status-title'),
  uploadStatusMessage: document.getElementById('upload-status-message'),
  uploadThumbnailStatus: document.getElementById('upload-thumbnail-status'),
  uploadPill: document.getElementById('upload-pill'),
  filtersForm: document.getElementById('filters-form'),
  pdbIdQueryInput: document.getElementById('pdb-id-query'),
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

function textInputValue(input) {
  const value = input ? String(input.value).trim() : '';
  if (value === '') return null;
  return value;
}

function normalizeSortBy(value) {
  const key = String(value || '').trim();
  return SORT_BY_ALIASES[key] || key;
}

function buildSearchParams() {
  const params = new URLSearchParams();
  const pdbIdQuery = textInputValue(elements.pdbIdQueryInput);
  if (pdbIdQuery !== null) params.set('pdb_id_query', pdbIdQuery);

  for (const field of elements.filtersForm.elements) {
    if (!field || !field.name) continue;
    if (field.name === 'sort_by' || field.name === 'sort_dir' || field.name === 'limit') {
      continue;
    }
    if ((field.type === 'checkbox' || field.type === 'radio') && !field.checked) {
      continue;
    }
    const value = String(field.value || '').trim();
    if (value === '') continue;
    params.set(field.name, value);
  }

  params.set('sort_by', normalizeSortBy(elements.sortBySelect.value));
  params.set('sort_dir', elements.sortDirSelect.value);
  params.set('limit', String(state.currentLimit));
  params.set('offset', String(state.currentOffset));
  return params;
}

async function fetchJson(url, options = {}) {
  const headers = new Headers(options.headers || {});
  if (!headers.has('Accept')) headers.set('Accept', 'application/json');
  const response = await fetch(url, { ...options, headers });
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
  state.health = health;
  elements.dbPill.textContent = health.db_exists ? health.db.split('/').slice(-2).join('/') : 'Missing DB';
  elements.uploadPill.textContent = health.upload_enabled ? 'Enabled' : 'Disabled';
  configureUploadForm(health);
}

function configureUploadForm(health) {
  const policies = Array.isArray(health.assembly_policies) && health.assembly_policies.length > 0
    ? health.assembly_policies
    : DEFAULT_ASSEMBLY_POLICIES;
  const previousPolicy = elements.uploadAssemblyPolicySelect.value;
  elements.uploadAssemblyPolicySelect.innerHTML = '';
  for (const policy of policies) {
    const option = document.createElement('option');
    option.value = String(policy);
    option.textContent = String(policy);
    elements.uploadAssemblyPolicySelect.append(option);
  }
  if (policies.includes(previousPolicy)) {
    elements.uploadAssemblyPolicySelect.value = previousPolicy;
  } else if (policies.includes('biological')) {
    elements.uploadAssemblyPolicySelect.value = 'biological';
  }

  elements.uploadResolutionInput.value = elements.uploadResolutionInput.value || DEFAULT_UPLOAD_RESOLUTION;
  elements.uploadForm.hidden = !health.upload_enabled;
  elements.uploadDisabledState.hidden = health.upload_enabled;
}

function clearUploadPoll() {
  if (state.uploadPollId !== null) {
    window.clearInterval(state.uploadPollId);
    state.uploadPollId = null;
  }
}

function setUploadBusy(isBusy) {
  for (const element of elements.uploadForm.elements) {
    element.disabled = isBusy;
  }
  elements.uploadResetButton.disabled = isBusy;
}

function formatBytes(bytes) {
  const numericBytes = Number(bytes);
  if (!Number.isFinite(numericBytes) || numericBytes <= 0) return '0 B';
  const units = ['B', 'KB', 'MB', 'GB'];
  let value = numericBytes;
  let unitIndex = 0;
  while (value >= 1024 && unitIndex < units.length - 1) {
    value /= 1024;
    unitIndex += 1;
  }
  return `${value.toFixed(value >= 10 || unitIndex === 0 ? 0 : 1)} ${units[unitIndex]}`;
}

function uploadStatusLabel(status) {
  const labels = {
    queued: 'Queued',
    running: 'Running analysis',
    indexing: 'Indexing result',
    done: 'Analysis complete',
    error: 'Analysis failed',
  };
  return labels[status] || 'Analyzing';
}

function thumbnailStatusLabel(status) {
  const labels = {
    pending: 'Thumbnail pending',
    rendering: 'Rendering thumbnail',
    done: 'Thumbnail ready',
    failed: 'Thumbnail failed',
    skipped: 'Thumbnail skipped',
  };
  return labels[status] || '';
}

function updateUploadProgressSteps(status) {
  const order = ['queued', 'running', 'indexing', 'done'];
  const activeIndex = status === 'error' ? -1 : Math.max(0, order.indexOf(status));
  for (const step of elements.uploadStatusPanel.querySelectorAll('[data-upload-step]')) {
    const stepIndex = order.indexOf(step.dataset.uploadStep);
    step.classList.toggle('is-active', stepIndex === activeIndex);
    step.classList.toggle('is-complete', activeIndex >= 0 && stepIndex < activeIndex);
    step.classList.toggle('is-error', status === 'error');
  }
}

function showUploadStatus(payload) {
  elements.uploadStatusPanel.hidden = false;
  const status = payload.status || 'queued';
  elements.uploadStatusTitle.textContent = uploadStatusLabel(status);
  elements.uploadStatusMessage.textContent = payload.error || payload.message || 'Waiting for the worker.';
  elements.uploadThumbnailStatus.textContent = thumbnailStatusLabel(payload.thumbnail_status);
  elements.uploadStatusPanel.classList.toggle('is-error', status === 'error');
  updateUploadProgressSteps(status);
}

async function pollUploadJob(jobId) {
  const payload = await fetchJson(`/api/analyze/${encodeURIComponent(jobId)}`);
  showUploadStatus(payload);

  if (payload.status === 'done' && payload.structure_id) {
    await loadRuns();
    await search();
    if (state.uploadRoutedJobId !== jobId) {
      state.uploadRoutedJobId = jobId;
      navigateToHit(payload.structure_id);
    }
    if (payload.thumbnail_status === 'done' || payload.thumbnail_status === 'failed' || payload.thumbnail_status === 'skipped') {
      clearUploadPoll();
      setUploadBusy(false);
    }
    return payload;
  }

  if (payload.status === 'error') {
    clearUploadPoll();
    setUploadBusy(false);
  }
  return payload;
}

async function submitUpload(event) {
  event.preventDefault();
  clearUploadPoll();
  state.uploadRoutedJobId = null;

  const file = elements.uploadFileInput.files[0];
  if (!file) {
    showUploadStatus({
      status: 'error',
      message: 'Choose a structure file before running analysis.',
      thumbnail_status: 'skipped',
    });
    return;
  }

  const maxUploadBytes = Number(state.health?.max_upload_bytes);
  if (Number.isFinite(maxUploadBytes) && maxUploadBytes > 0 && file.size > maxUploadBytes) {
    showUploadStatus({
      status: 'error',
      message: `Selected file is larger than ${formatBytes(maxUploadBytes)}.`,
      thumbnail_status: 'skipped',
    });
    return;
  }

  const formData = new FormData(elements.uploadForm);
  setUploadBusy(true);
  showUploadStatus({
    status: 'queued',
    message: 'Submitting upload.',
    thumbnail_status: 'pending',
  });

  try {
    const submitted = await fetchJson('/api/analyze', {
      method: 'POST',
      body: formData,
    });
    showUploadStatus({
      ...submitted,
      message: 'Queued for analysis.',
      thumbnail_status: 'pending',
    });
    const firstPoll = await pollUploadJob(submitted.job_id);
    const terminalThumbnail = ['done', 'failed', 'skipped'].includes(firstPoll.thumbnail_status);
    if (state.uploadPollId === null && firstPoll.status !== 'error' && !(firstPoll.status === 'done' && terminalThumbnail)) {
      state.uploadPollId = window.setInterval(() => {
        pollUploadJob(submitted.job_id).catch((error) => {
          showUploadStatus({
            status: 'error',
            message: `Unable to poll analysis job: ${error.message}`,
            thumbnail_status: 'skipped',
          });
          clearUploadPoll();
          setUploadBusy(false);
        });
      }, UPLOAD_POLL_INTERVAL_MS);
    }
  } catch (error) {
    showUploadStatus({
      status: 'error',
      message: error.message,
      thumbnail_status: 'skipped',
    });
    setUploadBusy(false);
  }
}

function resetUploadForm() {
  clearUploadPoll();
  state.uploadRoutedJobId = null;
  elements.uploadForm.reset();
  elements.uploadResolutionInput.value = DEFAULT_UPLOAD_RESOLUTION;
  configureUploadForm(state.health || {});
  elements.uploadStatusPanel.hidden = true;
  elements.uploadStatusPanel.classList.remove('is-error');
  setUploadBusy(false);
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
    openButton.addEventListener('click', () => navigateToHit(row.structure_id));
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

function appendTableCell(row, value) {
  const cell = document.createElement('td');
  cell.textContent = value;
  row.append(cell);
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
  const sortedVolumes = [...detail.volumes]
    .filter((volume) => String(volume.kind || '').toLowerCase() !== 'hub')
    .sort((a, b) => (b.volume_a3 ?? 0) - (a.volume_a3 ?? 0));
  for (const volume of sortedVolumes) {
    const row = document.createElement('tr');
    appendTableCell(row, volume.kind);
    appendTableCell(row, prettyNumber(volume.rank_in_kind, 0));
    appendTableCell(row, prettyNumber(volume.volume_a3, 0));
    appendTableCell(row, prettyNumber(volume.length_a));
    appendTableCell(row, prettyNumber(volume.min_diameter_a));
    appendTableCell(row, prettyNumber(volume.max_diameter_a));
    appendTableCell(row, prettyNumber(volume.cross_section_circularity, 2));
    appendTableCell(row, prettyNumber(volume.cross_section_uniformity, 2));
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
  await applyVolumeStyle(viewer, viewerData.volume_surface);
}

function getVolumeSurfaceTypeParams(volumeSurface) {
  const radiusOffset = Number(volumeSurface?.radius_offset);
  const smoothness = Number(volumeSurface?.smoothness);
  const quality = typeof volumeSurface?.quality === 'string' ? volumeSurface.quality : 'custom';

  return {
    quality,
    radiusOffset: Number.isFinite(radiusOffset) ? radiusOffset : 0,
    smoothness: Number.isFinite(smoothness) ? smoothness : 1.5,
  };
}

async function applyVolumeStyle(viewer, volumeSurface) {
  try {
    const plugin = viewer.plugin;
    if (!plugin || !plugin.managers || !plugin.managers.structure) return;

    const structures = plugin.managers.structure.hierarchy.current.structures;
    if (!structures || structures.length === 0) return;

    const structRef = structures[0];

    // Remove all default components and their representations
    const build = plugin.state.data.build();
    for (const comp of structRef.components) {
      build.delete(comp.cell.transform.ref);
    }
    await build.commit();

    // Add protein with cartoon representation in green/grey shades.
    // Use an explicit MolScript "not volume" selection instead of the
    // built-in 'polymer' preset, because Mol* cannot infer polymer
    // status without _entity_poly / _struct_asym tables (which biotite
    // does not write).
    const notVolumeExpr = NON_PROTEIN_COMP_IDS
      .map(id => `(= atom.label_comp_id ${id})`)
      .join(' ');
    const polymer = await plugin.builders.structure.tryCreateComponent(
      structRef.cell,
      {
        type: {
          name: 'script',
          params: {
            language: 'mol-script',
            expression: `(sel.atom.atom-groups :residue-test (not (or ${notVolumeExpr})))`,
          },
        },
        nullIfEmpty: true,
        label: 'Protein',
      },
      'protein',
    );
    if (polymer) {
      await plugin.builders.structure.representation.addRepresentation(
        polymer, {
          type: 'cartoon',
          color: 'chain-id',
          colorParams: {
            palette: {
              name: 'colors',
              params: {
                list: { kind: 'set', colors: PROTEIN_PALETTE },
              },
            },
          },
        },
      );
    }

    // Add per-type volume components with darker overlays for mouth shells.
    const volumeSurfaceTypeParams = getVolumeSurfaceTypeParams(volumeSurface);
    for (const [compId, style] of Object.entries(VOLUME_STYLES)) {
      const componentSpecs = [
        {
          key: `volume-${compId.toLowerCase()}-core`,
          selection: volumeScript(compId, 'core'),
          color: style.color,
        },
        {
          key: `volume-${compId.toLowerCase()}-mouth`,
          selection: volumeScript(compId, 'mouth'),
          color: style.mouthColor,
        },
      ];

      try {
        for (const spec of componentSpecs) {
          const comp = await plugin.builders.structure.tryCreateComponent(
            structRef.cell,
            spec.selection,
            spec.key,
          );
          if (comp) {
            await plugin.builders.structure.representation.addRepresentation(
              comp, {
                type: 'gaussian-surface',
                typeParams: volumeSurfaceTypeParams,
                color: 'uniform',
                colorParams: { value: spec.color },
              },
            );
          }
        }
      } catch (err) {
        console.warn(`Failed to create ${compId} volume component:`, err);
      }
    }
  } catch (error) {
    console.warn('Volume surface styling failed:', error);
  }
}

function getRoute() {
  const hash = window.location.hash || '#/';
  const hitMatch = hash.match(/^#\/hit\/(\d+)$/);
  if (hitMatch) {
    return { view: 'detail', structureId: Number(hitMatch[1]) };
  }
  if (hash === '#/gallery') {
    return { view: 'gallery', structureId: null };
  }
  return { view: 'landing', structureId: null };
}

function setCurrentView(view) {
  state.currentView = view;
  elements.landingView.hidden = view !== 'landing';
  elements.galleryView.hidden = view === 'landing';
  elements.homeViewButton.classList.toggle('is-active', view === 'landing');
  elements.galleryViewButton.classList.toggle('is-active', view !== 'landing');
}

function navigateTo(hash) {
  if (window.location.hash === hash) {
    renderRoute().catch((error) => {
      console.warn('Route update failed:', error);
    });
    return;
  }
  window.location.hash = hash;
}

function navigateToHit(structureId) {
  navigateTo(`#/hit/${structureId}`);
}

async function renderRoute() {
  const route = getRoute();
  if (route.view === 'landing') {
    setCurrentView('landing');
    if (elements.detailDialog.open) closeDetail({ updateHash: false });
    return;
  }

  setCurrentView('gallery');
  if (route.view === 'detail' && route.structureId) {
    await openDetail(route.structureId, { updateHash: false });
  } else if (elements.detailDialog.open) {
    closeDetail({ updateHash: false });
  }
}

async function openDetail(structureId, options = {}) {
  if (options.updateHash !== false && window.location.hash !== `#/hit/${structureId}`) {
    navigateToHit(structureId);
    return;
  }

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

function closeDetail(options = {}) {
  if (elements.detailDialog.open) {
    if (options.updateHash === false) {
      state.suppressDetailCloseRoute = true;
    }
    elements.detailDialog.close();
  }
  if (options.updateHash !== false && getRoute().view === 'detail') {
    navigateTo('#/gallery');
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

function loadStoredString(storageKey) {
  try {
    const value = localStorage.getItem(storageKey);
    return value === null ? '' : String(value);
  } catch {
    return '';
  }
}

function saveStoredString(storageKey, value) {
  if (!value) {
    localStorage.removeItem(storageKey);
    return;
  }
  localStorage.setItem(storageKey, String(value));
}

function createPresetOption(value, label) {
  const opt = document.createElement('option');
  opt.value = value;
  opt.textContent = label;
  return opt;
}

function hasSelectOption(selectEl, value) {
  return [...selectEl.options].some((option) => option.value === value);
}

function populatePresetSelect(selectEl, options) {
  const previousValue = selectEl.value;
  selectEl.innerHTML = '';
  for (const option of options) {
    selectEl.append(createPresetOption(option.value, option.label));
  }
  if (hasSelectOption(selectEl, previousValue)) {
    selectEl.value = previousValue;
  }
}

function getBuiltinFilterPresetValue(name) {
  return `${BUILTIN_FILTER_PRESET_PREFIX}${name}`;
}

function isBuiltinFilterPresetValue(value) {
  return String(value || '').startsWith(BUILTIN_FILTER_PRESET_PREFIX);
}

function getBuiltinFilterPresetName(value) {
  if (!isBuiltinFilterPresetValue(value)) return '';
  return String(value).slice(BUILTIN_FILTER_PRESET_PREFIX.length);
}

function isReservedFilterPresetName(name) {
  return name === FILTER_PRESET_NONE_LABEL || Object.hasOwn(BUILTIN_FILTER_PRESETS, name);
}

function loadUserFilterPresetsMap() {
  const map = loadPresetsMap(FILTER_PRESETS_KEY);
  const filteredMap = {};
  for (const [name, preset] of Object.entries(map)) {
    if (!isReservedFilterPresetName(name)) {
      filteredMap[name] = preset;
    }
  }
  return filteredMap;
}

function getFilterPresetData(value) {
  if (!value) return {};
  if (isBuiltinFilterPresetValue(value)) {
    const builtinName = getBuiltinFilterPresetName(value);
    return BUILTIN_FILTER_PRESETS[builtinName] || null;
  }
  const map = loadUserFilterPresetsMap();
  return map[value] || null;
}

function populateFilterPresetSelect() {
  const options = [{ value: '', label: FILTER_PRESET_NONE_LABEL }];
  for (const name of Object.keys(BUILTIN_FILTER_PRESETS)) {
    options.push({ value: getBuiltinFilterPresetValue(name), label: name });
  }
  for (const name of Object.keys(loadUserFilterPresetsMap()).sort()) {
    options.push({ value: name, label: name });
  }
  populatePresetSelect(elements.filterPresetSelect, options);
}

function getBuiltinDisplayPresetValue(name) {
  return `${BUILTIN_DISPLAY_PRESET_PREFIX}${name}`;
}

function isBuiltinDisplayPresetValue(value) {
  return String(value || '').startsWith(BUILTIN_DISPLAY_PRESET_PREFIX);
}

function getBuiltinDisplayPresetName(value) {
  if (!isBuiltinDisplayPresetValue(value)) return '';
  return String(value).slice(BUILTIN_DISPLAY_PRESET_PREFIX.length);
}

function isReservedDisplayPresetName(name) {
  return name === DISPLAY_PRESET_DEFAULT_LABEL || Object.hasOwn(BUILTIN_DISPLAY_PRESETS, name);
}

function loadUserDisplayPresetsMap() {
  const map = loadPresetsMap(DISPLAY_PRESETS_KEY);
  const filteredMap = {};
  for (const [name, preset] of Object.entries(map)) {
    if (!isReservedDisplayPresetName(name)) {
      filteredMap[name] = preset;
    }
  }
  return filteredMap;
}

function getDisplayPresetData(value) {
  if (!value) return DEFAULT_DISPLAY_METRICS;
  if (isBuiltinDisplayPresetValue(value)) {
    const builtinName = getBuiltinDisplayPresetName(value);
    return BUILTIN_DISPLAY_PRESETS[builtinName] || null;
  }
  const map = loadUserDisplayPresetsMap();
  return map[value] || null;
}

function populateDisplayPresetSelect() {
  const options = [{ value: '', label: DISPLAY_PRESET_DEFAULT_LABEL }];
  for (const name of Object.keys(BUILTIN_DISPLAY_PRESETS)) {
    options.push({ value: getBuiltinDisplayPresetValue(name), label: name });
  }
  for (const name of Object.keys(loadUserDisplayPresetsMap()).sort()) {
    options.push({ value: name, label: name });
  }
  populatePresetSelect(elements.displayPresetSelect, options);
}

function getActiveFilterPresetName() {
  return loadStoredString(ACTIVE_FILTER_PRESET_KEY).trim();
}

function setActiveFilterPresetName(name) {
  saveStoredString(ACTIVE_FILTER_PRESET_KEY, String(name || '').trim());
}

function captureFilterState() {
  const data = {};
  for (const el of elements.filtersForm.elements) {
    if (el.name) {
      data[el.name] = el.name === 'sort_by' ? normalizeSortBy(el.value) : el.value;
    }
  }
  data.pdb_id_query = elements.pdbIdQueryInput.value;
  return data;
}

function clearFilterState() {
  elements.filtersForm.reset();
  for (const el of elements.filtersForm.elements) {
    if (!el || !el.name) continue;
    if (el.name === 'sort_by') {
      el.value = normalizeSortBy(el.value);
      continue;
    }
    if (NON_FILTER_FORM_FIELD_NAMES.has(el.name)) continue;
    if (el.type === 'checkbox' || el.type === 'radio') {
      el.checked = false;
      continue;
    }
    el.value = '';
  }
  elements.pdbIdQueryInput.value = '';
}

function applyFilterState(data) {
  // Reset sort/page controls to their form defaults, then explicitly blank
  // actual filters so browser-restored state cannot survive a "None" load.
  clearFilterState();
  for (const el of elements.filtersForm.elements) {
    if (el.name && data[el.name] !== undefined) {
      el.value = el.name === 'sort_by' ? normalizeSortBy(data[el.name]) : data[el.name];
    }
  }
  elements.pdbIdQueryInput.value = data.pdb_id_query ?? '';
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

function restoreActiveFilterPreset() {
  const name = getActiveFilterPresetName();
  if (!name) return false;

  const preset = getFilterPresetData(name);
  if (!preset) {
    setActiveFilterPresetName('');
    return false;
  }

  if (hasSelectOption(elements.filterPresetSelect, name)) {
    elements.filterPresetSelect.value = name;
  }
  applyFilterState(preset);
  return true;
}

function promptPresetName(action, suggestedName = '') {
  const name = prompt(`Enter a name to ${action}:`, suggestedName);
  if (!name || !name.trim()) return null;
  return name.trim();
}

function wirePresets() {
  // Filter presets
  populateFilterPresetSelect();
  const activeFilterPresetName = getActiveFilterPresetName();
  if (
    activeFilterPresetName &&
    hasSelectOption(elements.filterPresetSelect, activeFilterPresetName)
  ) {
    elements.filterPresetSelect.value = activeFilterPresetName;
  }

  elements.filterPresetSave.addEventListener('click', () => {
    const selectedValue = elements.filterPresetSelect.value;
    const name = promptPresetName(
      'save this filter preset',
      isBuiltinFilterPresetValue(selectedValue) ? '' : selectedValue,
    );
    if (!name) return;
    if (isReservedFilterPresetName(name)) {
      alert(`"${name}" is reserved for built-in presets.`);
      return;
    }
    const map = loadUserFilterPresetsMap();
    if (map[name] && !confirm(`Overwrite existing preset "${name}"?`)) return;
    map[name] = captureFilterState();
    savePresetsMap(FILTER_PRESETS_KEY, map);
    setActiveFilterPresetName(name);
    populateFilterPresetSelect();
    elements.filterPresetSelect.value = name;
  });

  elements.filterPresetLoad.addEventListener('click', () => {
    const name = elements.filterPresetSelect.value;
    const preset = getFilterPresetData(name);
    if (!preset) return;
    setActiveFilterPresetName(name);
    applyFilterState(preset);
    search();
  });

  elements.filterPresetDelete.addEventListener('click', () => {
    const name = elements.filterPresetSelect.value;
    if (!name || isBuiltinFilterPresetValue(name)) return;
    if (!confirm(`Delete preset "${name}"?`)) return;
    const map = loadUserFilterPresetsMap();
    delete map[name];
    savePresetsMap(FILTER_PRESETS_KEY, map);
    if (getActiveFilterPresetName() === name) {
      setActiveFilterPresetName('');
    }
    populateFilterPresetSelect();
  });

  // Display property presets
  populateDisplayPresetSelect();

  elements.displayPresetSave.addEventListener('click', () => {
    const selectedValue = elements.displayPresetSelect.value;
    const name = promptPresetName(
      'save this display preset',
      isBuiltinDisplayPresetValue(selectedValue) ? '' : selectedValue,
    );
    if (!name) return;
    if (isReservedDisplayPresetName(name)) {
      alert(`"${name}" is reserved for built-in presets.`);
      return;
    }
    const map = loadUserDisplayPresetsMap();
    if (map[name] && !confirm(`Overwrite existing preset "${name}"?`)) return;
    map[name] = captureDisplayState();
    savePresetsMap(DISPLAY_PRESETS_KEY, map);
    populateDisplayPresetSelect();
    elements.displayPresetSelect.value = name;
  });

  elements.displayPresetLoad.addEventListener('click', () => {
    const name = elements.displayPresetSelect.value;
    const preset = getDisplayPresetData(name);
    if (!preset) return;
    applyDisplayState(preset);
    search();
  });

  elements.displayPresetDelete.addEventListener('click', () => {
    const name = elements.displayPresetSelect.value;
    if (!name || isBuiltinDisplayPresetValue(name)) return;
    if (!confirm(`Delete preset "${name}"?`)) return;
    const map = loadUserDisplayPresetsMap();
    delete map[name];
    savePresetsMap(DISPLAY_PRESETS_KEY, map);
    populateDisplayPresetSelect();
  });
}

function wireEvents() {
  elements.homeViewButton.addEventListener('click', () => navigateTo('#/'));
  elements.galleryViewButton.addEventListener('click', () => navigateTo('#/gallery'));
  elements.browseGalleryButton.addEventListener('click', () => navigateTo('#/gallery'));
  elements.uploadFocusButton.addEventListener('click', () => {
    navigateTo('#/');
    window.setTimeout(() => elements.uploadFileInput.focus(), 0);
  });
  elements.uploadForm.addEventListener('submit', submitUpload);
  elements.uploadResetButton.addEventListener('click', resetUploadForm);
  window.addEventListener('hashchange', () => {
    renderRoute().catch((error) => {
      console.warn('Route update failed:', error);
    });
  });

  elements.filtersForm.addEventListener('submit', async (event) => {
    event.preventDefault();
    state.currentOffset = 0;
    state.currentLimit = Number(elements.limitSelect.value);
    await search();
  });

  elements.pdbIdQueryInput.addEventListener('keydown', async (event) => {
    if (event.key !== 'Enter') return;
    event.preventDefault();
    state.currentOffset = 0;
    state.currentLimit = Number(elements.limitSelect.value);
    await search();
  });

  elements.resetButton.addEventListener('click', async (event) => {
    event.preventDefault();

    const restoredActivePreset = restoreActiveFilterPreset();
    if (!restoredActivePreset) {
      clearFilterState();
    }

    state.currentLimit = Number(elements.limitSelect.value);
    state.currentOffset = 0;
    await search();
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
  elements.detailDialog.addEventListener('close', () => {
    if (state.suppressDetailCloseRoute) {
      state.suppressDetailCloseRoute = false;
      return;
    }
    if (getRoute().view === 'detail') {
      navigateTo('#/gallery');
    }
  });
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
  restoreActiveFilterPreset();
  await search();
  await renderRoute();
}

window.addEventListener('DOMContentLoaded', () => {
  init().catch((error) => {
    elements.loadingState.hidden = true;
    elements.emptyState.hidden = false;
    elements.emptyState.textContent = `Startup failed: ${error.message}`;
  });
});
