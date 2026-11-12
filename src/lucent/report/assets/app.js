/* LUCENT offline evidence explorer. No network requests or inferential calculations. */
(function (root) {
  'use strict';
  const finite = value => typeof value === 'number' && Number.isFinite(value);
  const probability = value => finite(value) && value >= 0 && value <= 1;
  const resultFor = (unit, version) => unit && unit.versions ? unit.versions[version] || null : null;
  const unavailableStatus = status => /degenerate|missing[_ -]?context|no[_ -]?structural[_ -]?edges|missing[_ -]?state|unavailable|ineligible|insufficient|failed|error|unsupported|not[_ -]?evaluable|\u4e0d\u53ef\u7528|\u4e0d\u8db3/i.test(String(status || ''));
  const evaluable = result => !!result && result.evaluable !== false && !unavailableStatus(result.status) && probability(result.p);
  function adjustment(version) {
    const method = String(version && version.adjustment_method || '').trim();
    if (/holm/i.test(method)) return { kind: 'holm', label: 'Holm adjusted p', highlight: 'Holm highlight', eligible: true };
    if (/^(BH|BY)$/i.test(method) || /benjamini/i.test(method)) {
      const name = /BY|yekutieli/i.test(method) ? 'BY' : 'BH';
      return { kind: 'fdr', label: name + ' FDR q', highlight: name + ' FDR highlight', eligible: true };
    }
    return { kind: 'none', label: method ? method + ' Adjusted value' : 'Adjusted value', highlight: 'Adjusted-result highlight', eligible: false };
  }
  function missingLabel(result) {
    if (!result) return 'No result in this version';
    const status = String(result.status || '');
    if (unavailableStatus(status) || result.evaluable === false && !/not[_ -]?computed|not[_ -]?run|pending|planned|\u672a\u8ba1\u7b97|\u672a\u8fd0\u884c/i.test(status)) return 'Not evaluable';
    return 'Not computed';
  }
  function displayStatus(status, record) {
    const raw = String(status || '');
    if (/MISSING_CONTEXT/i.test(raw)) return 'Missing context';
    if (/NO_STRUCTURAL_EDGES/i.test(raw)) return 'No spatial edges';
    if (/DEGENERATE/i.test(raw)) return 'Degenerate expression';
    if (/RESPONSE_CONSTANT/i.test(raw)) return 'Constant response';
    if (unavailableStatus(raw) || record && record.evaluable === false) return missingLabel(record || {status:raw});
    if (/not[_ -]?computed|not[_ -]?run|pending|planned/i.test(raw)) return 'Not computed';
    if (/historical|legacy/i.test(raw)) return 'Historical';
    if (/sensitivity/i.test(raw)) return 'Sensitivity analysis';
    if (/research|exploratory/i.test(raw)) return 'Exploratory';
    if (/available|evaluable|computed|conditional_support/i.test(raw)) return 'Available';
    return raw ? 'Recorded' : 'Status unspecified';
  }
  function formatNumber(value, kind) {
    if (!finite(value)) return '—';
    if (value === 0) return '0';
    if (kind === 'p' || kind === 'q') return Math.abs(value) < .001 ? value.toExponential(2) : value.toPrecision(3).replace(/0+$/, '').replace(/\.$/, '');
    return Math.abs(value) >= 10000 || Math.abs(value) < .0001 ? value.toExponential(2) : Number(value.toPrecision(4)).toString();
  }
  const exact = value => finite(value) ? String(value) : 'null';
  function evidencePass(result, mode, version) {
    if (mode === 'all') return true;
    if (mode === 'evaluable') return evaluable(result);
    if (mode === 'nominal') return evaluable(result) && result.p < .05;
    if (mode === 'adjusted') return evaluable(result) && adjustment(version).eligible && probability(result.q) && result.q < .05;
    return false;
  }
  function geneticSearchIndex(links) {
    const index = new Map();
    for (const link of links || []) index.set(link.channel_id, (index.get(link.channel_id) || '') + ' ' + [link.locus_id, link.variant_id, link.rsid, link.gene_id, link.gene_symbol].filter(Boolean).join(' ').toLocaleLowerCase());
    return index;
  }
  function filterUnits(units, filters, version, geneticLinks) {
    const needle = String(filters.search || '').trim().toLocaleLowerCase();
    const genetics = geneticLinks instanceof Map ? geneticLinks : geneticSearchIndex(geneticLinks);
    return units.filter(unit => (!filters.context || unit.context === filters.context) && (!filters.sender || unit.sender_state === filters.sender) && (!filters.receiver || unit.receiver_state === filters.receiver) && (!needle || ([unit.unit_id, unit.channel_id, unit.lr_pair, unit.context, unit.sender_state, unit.receiver_state].concat(unit.ligand_genes || [], unit.receptor_genes || []).join(' ').toLocaleLowerCase() + ' ' + (genetics.get(unit.channel_id) || '')).includes(needle)) && evidencePass(resultFor(unit, version.id), filters.evidence || 'all', version));
  }
  function sortUnits(units, versionId, key, ascending) {
    return units.slice().sort((a, b) => {
      const av = (resultFor(a, versionId) || {})[key];
      const bv = (resultFor(b, versionId) || {})[key];
      if (!finite(av) && !finite(bv)) return String(a.unit_id).localeCompare(String(b.unit_id));
      if (!finite(av)) return 1;
      if (!finite(bv)) return -1;
      return (ascending ? av - bv : bv - av) || String(a.unit_id).localeCompare(String(b.unit_id));
    });
  }
  const ramp = [[59,130,156],[88,168,172],[145,184,160],[213,194,127],[220,142,85],[195,82,68]];
  function pColor(p) {
    if (!probability(p)) return '#a6ada8';
    const scaled = Math.min(6, p === 0 ? 6 : -Math.log10(p)) / 6 * (ramp.length - 1);
    const index = Math.min(ramp.length - 2, Math.floor(scaled));
    const f = scaled - index;
    return 'rgb(' + ramp[index].map((v, j) => Math.round(v + f * (ramp[index + 1][j] - v))).join(',') + ')';
  }
  function orientedEdges(section, unit) {
    if (!section || !unit) return [];
    if (section.context != null && section.context !== unit.context) return [];
    const output = [], seen = new Set();
    for (const edge of section.edges || []) {
      const a = section.spots[edge[0]], b = section.spots[edge[1]];
      if (!a || !b) continue;
      let from, to;
      if (a.state === unit.sender_state && b.state === unit.receiver_state) { from = edge[0]; to = edge[1]; }
      else continue;
      const key = from + ':' + to;
      if (!seen.has(key)) { output.push([from, to]); seen.add(key); }
    }
    return output;
  }
  function tsvCell(value) {
    if (value === null || value === undefined) return 'null';
    let text = String(value).replace(/[\t\r\n]+/g, ' ');
    if (/^[=+@]/.test(text) || /^-[^\d.]/.test(text)) text = "'" + text;
    return text;
  }
  function unitsTSV(units, version) {
    const columns = ['unit_id','channel_id','lr_pair','context','sender_state','receiver_state','ligand_genes','receptor_genes','version','score','p','q','p_label','q_label','adjustment_method','family_size','status','evaluable'];
    const rows = units.map(unit => { const r = resultFor(unit, version.id) || {}; return [unit.unit_id,unit.channel_id,unit.lr_pair,unit.context,unit.sender_state,unit.receiver_state,(unit.ligand_genes || []).join(';'),(unit.receptor_genes || []).join(';'),version.id,r.score,r.p,r.q,version.p_label || 'p',adjustment(version).kind === 'holm' ? 'Holm adjusted p' : version.q_label || adjustment(version).label,version.adjustment_method,version.family_size,r.status,r.evaluable].map(tsvCell).join('\t'); });
    return columns.join('\t') + '\r\n' + rows.join('\r\n') + '\r\n';
  }
  function fitScaleForBounds(bounds, width, height) {
    return Math.min(Math.max(1, width - 88) / Math.max(bounds.xmax - bounds.xmin, 1), Math.max(1, height - 88) / Math.max(bounds.ymax - bounds.ymin, 1));
  }
  function resizedView(view, nextFitScale) {
    const zoom = finite(view.scale / view.fitScale) ? Math.max(.15, Math.min(100, view.scale / view.fitScale)) : 1;
    return { ...view, fitScale: nextFitScale, scale: nextFitScale * zoom };
  }
  function pointKindInfo(section, metadata, unit) {
    const description = String(section && section.point_kind || metadata && metadata.point_kind || unit && unit.point_kind || 'Spatial spots');
    const bins = /bin/i.test(description);
    return { description, singular: bins ? 'Spatial bin' : 'Spot', plural: bins ? 'spatial bins' : 'spots', states: bins ? 'Spatial bin states' : 'Cell / spot states' };
  }
  function savedDiagnostic(result, key) {
    if (!result) return null;
    if (Object.prototype.hasOwnProperty.call(result, key)) return result[key];
    return result.source_fields && Object.prototype.hasOwnProperty.call(result.source_fields, key) ? result.source_fields[key] : null;
  }
  const sectionEvidence = metadata => /^(single_embryo_section|paired_embryo_sections)$/.test(metadata.analysis_scope || '');
  const replicateLabel = metadata => metadata.replicate_label || 'Donor';
  const logic = { finite, probability, resultFor, evaluable, unavailableStatus, adjustment, missingLabel, displayStatus, formatNumber, evidencePass, geneticSearchIndex, filterUnits, sortUnits, pColor, orientedEdges, tsvCell, unitsTSV, fitScaleForBounds, resizedView, pointKindInfo, savedDiagnostic, sectionEvidence, replicateLabel };
  root.LUCENT_LOGIC = logic;
  if (typeof module === 'object' && module.exports) module.exports = logic;
  if (typeof document === 'undefined') return;

  const $ = id => document.getElementById(id);
  function node(tag, className, text) { const n = document.createElement(tag); if (className) n.className = className; if (text !== undefined && text !== null) n.textContent = String(text); return n; }
  function setText(id, text) { $(id).textContent = text == null ? '' : String(text); }
  function empty(element) { element.replaceChildren(); }
  function option(value, label) { const o = node('option', null, label); o.value = value; return o; }
  function fillSelect(id, values, allLabel, chosen) { const el = $(id); empty(el); if (allLabel !== null) el.append(option('', allLabel)); for (const value of values) el.append(option(typeof value === 'object' ? value.id : value, typeof value === 'object' ? value.label : value)); el.value = chosen || ''; if (el.selectedIndex < 0 && el.options.length) el.selectedIndex = 0; return el.value; }
  const unique = values => [...new Set(values.filter(v => v !== null && v !== undefined && v !== ''))].sort((a,b) => String(a).localeCompare(String(b), undefined, { numeric: true }));
  const report = root.LUCENT_REPORT;
  if (!report || !report.metadata || !Array.isArray(report.units)) { $('global-error').hidden = false; setText('global-error', 'Report data could not be loaded. Keep index.html, assets, and data in their original relative locations, then reopen index.html.'); return; }
  report.genetic_links = report.genetic_links || [];
  report.donor_evidence = report.donor_evidence || [];
  report.sections = report.sections || [];
  const geneticSearch = geneticSearchIndex(report.genetic_links);
  const metadata = report.metadata;
  const versions = metadata.versions || [];
  if (!versions.length) { $('global-error').hidden = false; setText('global-error', 'Result-version metadata is missing from this report.'); return; }
  const state = {
    version: versions.find(v => v.id === metadata.default_version) || versions[0],
    filters: { context: '', sender: '', receiver: '', search: '', evidence: 'all' },
    selected: null, sectionMeta: null, section: null, image: null,
    filtered: [], page: 0, pageSize: 50, sort: 'p', ascending: true,
    mode: 'pan', gene: '', showEdges: true, highlight: false, showImage: true,
    view: null, bounds: null, dragging: null, edges: [], loadToken: 0, hovered: null
  };
  const canvas = $('spatial-canvas'), ctx = canvas.getContext('2d');
  let width = 500, height = 430, drawPending = false, toastTimer = 0;
  const sectionPromises = new Map();
  const colors = ['#73999a','#bc987b','#a2b082','#988aaf','#6d9dba','#bc91a0','#d3b471','#83afa0','#a6a49a','#b29eaa','#99aec2','#c3ad95'];
  const stateColors = new Map(unique(report.units.flatMap(u => [u.sender_state, u.receiver_state])).map((name,i) => [name, colors[i % colors.length]]));
  function stateColor(name) { if (!stateColors.has(name)) stateColors.set(name, colors[stateColors.size % colors.length]); return stateColors.get(name); }
  function currentResult() { return resultFor(state.selected, state.version.id); }
  function currentPointKind() { return pointKindInfo(state.sectionMeta && state.sectionMeta.point_kind ? state.sectionMeta : state.section, metadata, state.selected); }
  function correctionLabel() { return adjustment(state.version).label; }
  function toast(text) { clearTimeout(toastTimer); setText('toast', text); $('toast').hidden = false; toastTimer = setTimeout(() => { $('toast').hidden = true; }, 2600); }
  function candidateLabel(unit) { return unit.lr_pair || [(unit.ligand_genes || []).join('+'), (unit.receptor_genes || []).join('+')].filter(Boolean).join(' → ') || unit.channel_id; }

  function initialize() {
    document.title = metadata.title || 'LUCENT Spatial Evidence Report';
    setText('report-title', metadata.title || 'Genetically Supported Spatial Communication Evidence');
    setText('report-subtitle', [metadata.dataset, metadata.trait].filter(Boolean).join(' · '));
    setText('claim', metadata.claim || 'Explore upstream genetic support alongside spatial candidate evidence, with one statistical record per hypothesis.');
    setText('total-units', report.units.length.toLocaleString());
    setText('section-count', report.sections.length.toLocaleString());
    setText('tab-donors', sectionEvidence(metadata) ? 'Section evidence' : 'Donor evidence');
    setText('replicate-label', replicateLabel(metadata));
    fillSelect('version', versions.map(v => ({ id: v.id, label: v.label || v.id })), null, state.version.id);
    fillSelect('context', unique(report.units.map(u => u.context).concat(report.sections.map(s => s.context))), 'All contexts', '');
    fillSelect('sender', unique(report.units.map(u => u.sender_state)), 'All sender states', '');
    fillSelect('receiver', unique(report.units.map(u => u.receiver_state)), 'All receiver states', '');
    state.filters.evidence = report.units.some(u => evaluable(resultFor(u, state.version.id))) ? 'evaluable' : 'all';
    renderSources();
    bindEvents();
    updateVersion();
    updateSections(metadata.default_section_id);
    updateCandidates(false);
    if (typeof ResizeObserver === 'function') new ResizeObserver(resizeCanvas).observe($('canvas-wrap'));
    else root.addEventListener('resize', resizeCanvas);
    resizeCanvas();
  }
  function renderSources() {
    const container = $('sources'); empty(container);
    if (metadata.source_status) container.append(node('p', 'source-label', metadata.source_status));
    for (const source of metadata.sources || []) {
      const p = node('p'); p.append(node('strong', null, source.label || 'Source'), node('br'), node('span', null, source.path || 'Path not provided'));
      if (source.sha256) { p.append(node('br'), node('code', null, 'SHA256 ' + source.sha256)); }
      container.append(p);
    }
    for (const note of metadata.notes || []) container.append(node('p', null, note));
    if (!container.childNodes.length) container.append(node('p', null, 'No additional source notes were provided.'));
  }
  function updateVersion() {
    const v = state.version, adj = adjustment(v);
    setText('evaluable-units', report.units.filter(u => evaluable(resultFor(u, v.id))).length.toLocaleString());
    setText('adjusted-filter-label', adj.kind === 'holm' ? 'Holm adjusted p < 0.05' : adj.kind === 'fdr' ? adj.label + ' < 0.05' : 'Adjusted < 0.05');
    setText('highlight-label', adj.highlight);
    $('highlight-adjusted').disabled = !adj.eligible;
    if (!adj.eligible) { state.highlight = false; $('highlight-adjusted').checked = false; }
    setText('p-heading', (v.p_label || 'p') + ' ↕');
    setText('q-heading', correctionLabel() + ' ↕');
    setText('version-note', [displayStatus(v.status), v.adjustment_method ? 'Adjustment: ' + v.adjustment_method + (v.family_size != null ? ' · N = ' + v.family_size : '') : 'Adjustment: not provided'].join('\n'));
    $('version-note').style.whiteSpace = 'pre-line';
    setText('footer-version', 'Version ' + v.id + ' · ' + (v.adjustment_method || 'Adjustment unspecified') + ' · schema ' + (metadata.schema_version || '—'));
    document.querySelectorAll('input[name=evidence]').forEach(el => { el.checked = el.value === state.filters.evidence; });
  }
  function baseFilters() { return { ...state.filters, evidence: 'all' }; }
  function updateCandidates(keepSelection = true) {
    const base = filterUnits(report.units, baseFilters(), state.version, geneticSearch);
    for (const mode of ['all', 'evaluable', 'nominal', 'adjusted']) setText('count-' + mode, base.filter(u => evidencePass(resultFor(u, state.version.id), mode, state.version)).length.toLocaleString());
    state.filtered = sortUnits(filterUnits(report.units, state.filters, state.version, geneticSearch), state.version.id, state.sort, state.ascending);
    state.page = Math.min(state.page, Math.max(0, Math.ceil(state.filtered.length / state.pageSize) - 1));
    if (!keepSelection || !state.selected || !state.filtered.some(u => u.unit_id === state.selected.unit_id)) state.selected = state.filtered.find(u => u.unit_id === metadata.default_unit_id) || state.filtered[0] || null;
    setText('candidate-count', state.filtered.length.toLocaleString());
    setText('table-caption', 'Each row is one channel × context × state-pair hypothesis. Version: ' + (state.version.label || state.version.id));
    renderTable(); renderEvidence(); refreshEdges();
    if (state.selected && (!state.sectionMeta || state.sectionMeta.context !== state.selected.context)) chooseSectionForUnit(state.selected);
  }
  function renderTable() {
    const body = $('candidate-rows'); empty(body);
    const offset = state.page * state.pageSize;
    for (const unit of state.filtered.slice(offset, offset + state.pageSize)) {
      const r = resultFor(unit, state.version.id), row = node('tr', state.selected && state.selected.unit_id === unit.unit_id ? 'selected' : '');
      row.tabIndex = 0; row.setAttribute('aria-selected', String(!!state.selected && state.selected.unit_id === unit.unit_id)); row.dataset.unitId = unit.unit_id;
      const first = node('td'); first.append(node('strong', null, candidateLabel(unit)), node('span', 'subtext', unit.channel_id)); row.append(first, node('td', null, unit.sender_state + ' → ' + unit.receiver_state), node('td', null, unit.context));
      for (const key of ['score','p','q']) { const value = r && r[key], showValue = finite(value) && (key === 'score' || evaluable(r)); const td = node('td', 'numeric' + (showValue ? '' : ' missing'), showValue ? formatNumber(value, key) : missingLabel(r)); td.title = key + ' = ' + exact(value) + (!showValue && finite(value) ? ' · source placeholder' : ''); row.append(td); }
      const status = node('td'), dot = node('i', 'status-dot' + (evaluable(r) ? '' : missingLabel(r) === 'Not evaluable' ? ' unavailable' : ' pending'));
      status.append(dot, document.createTextNode(r ? displayStatus(r.status, r) : missingLabel(r))); row.append(status);
      row.addEventListener('click', () => selectUnit(unit));
      row.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); selectUnit(unit); } });
      body.append(row);
    }
    $('table-empty').hidden = state.filtered.length > 0;
    setText('table-range', state.filtered.length ? 'Showing ' + (offset + 1) + '–' + Math.min(offset + state.pageSize, state.filtered.length) + ' / ' + state.filtered.length : '0 candidates');
    setText('page-label', (state.page + 1) + ' / ' + Math.max(1, Math.ceil(state.filtered.length / state.pageSize)));
    $('page-prev').disabled = state.page === 0; $('page-next').disabled = offset + state.pageSize >= state.filtered.length;
    $('export-tsv').disabled = state.filtered.length === 0;
  }
  function selectUnit(unit) {
    state.selected = unit; renderTable(); renderEvidence();
    chooseSectionForUnit(unit); refreshEdges();
  }
  function availableSections() { return report.sections.filter(s => !state.filters.context || s.context === state.filters.context); }
  function updateSections(preferred) {
    const sections = availableSections();
    const chosen = sections.find(s => s.id === preferred);
    let donor = chosen ? chosen.donor_id : $('donor').value;
    donor = fillSelect('donor', unique(sections.map(s => s.donor_id)), 'All ' + replicateLabel(metadata).toLowerCase() + 's', donor);
    const subset = sections.filter(s => !donor || s.donor_id === donor);
    const sectionId = fillSelect('section', subset.map(s => ({ id: s.id, label: s.label || s.id })), subset.length ? null : 'No matching section', preferred || (state.sectionMeta && state.sectionMeta.id));
    const next = report.sections.find(s => s.id === sectionId) || null;
    if (next && state.sectionMeta && next.id === state.sectionMeta.id && state.section) return;
    loadSection(next);
  }
  function chooseSectionForUnit(unit) {
    if (!unit) return;
    const matching = availableSections().filter(s => s.context === unit.context);
    if (state.sectionMeta && matching.some(s => s.id === state.sectionMeta.id)) return;
    const next = matching.find(s => s.donor_id === $('donor').value) || matching[0];
    if (next) updateSections(next.id);
    else loadSection(null);
  }
  function sectionScript(meta) {
    root.LUCENT_SECTIONS = root.LUCENT_SECTIONS || {};
    if (root.LUCENT_SECTIONS[meta.id]) return Promise.resolve(root.LUCENT_SECTIONS[meta.id]);
    if (sectionPromises.has(meta.id)) return sectionPromises.get(meta.id);
    const promise = new Promise((resolve, reject) => { const script = document.createElement('script'); script.src = meta.file; script.async = true; script.onload = () => root.LUCENT_SECTIONS[meta.id] ? resolve(root.LUCENT_SECTIONS[meta.id]) : reject(new Error('The section script does not contain the expected data')); script.onerror = () => reject(new Error('Could not load the local section file: ' + meta.file)); document.body.append(script); });
    sectionPromises.set(meta.id, promise); return promise;
  }
  async function loadSection(meta) {
    const token = ++state.loadToken;
    state.sectionMeta = meta; state.section = null; state.image = null; state.view = null; state.bounds = null; state.edges = []; state.hovered = null;
    $('canvas-tooltip').hidden = true; $('canvas-status').hidden = false;
    if (!meta) { setText('canvas-status', 'No section matches the current selection. Candidate statistics remain available in the table.'); setText('section-title', 'No matching section'); setText('section-description', ''); setText('coordinate-label', 'Coordinates'); empty($('state-legend')); setText('plot-summary', ''); requestDraw(); return; }
    setText('canvas-status', 'Loading section…');
    setText('section-title', meta.label || meta.id);
    const pointKind = pointKindInfo(meta, metadata, state.selected);
    setText('section-description', [meta.context, replicateLabel(metadata) + ' ' + meta.donor_id, Number(meta.n_spots || 0).toLocaleString() + ' ' + pointKind.plural, Number(meta.n_edges || 0).toLocaleString() + ' source adjacency edges'].join(' · '));
    canvas.setAttribute('aria-label', pointKind.description + ' and directed candidate channel edges. Scroll to zoom; drag to pan or select a region.');
    const coordinate = meta.coordinate_system || {};
    setText('coordinate-label', (coordinate.unit || 'Source coordinates') + ' · y ' + (coordinate.y_direction === 'up' ? '↑' : '↓'));
    try {
      const section = await sectionScript(meta);
      if (token !== state.loadToken) return;
      state.section = section;
      for (const spot of section.spots || []) stateColor(spot.state);
      if (!section.spots || !section.spots.length) { setText('canvas-status', 'This section has no spots to display.'); return; }
      $('canvas-status').hidden = true;
      const genes = Object.keys(section.expression || {}).filter(gene => Array.isArray(section.expression[gene]) && section.expression[gene].length === section.spots.length);
      empty($('color-mode')); $('color-mode').append(option('state', currentPointKind().states));
      for (const gene of genes) $('color-mode').append(option('gene:' + gene, 'Expression · ' + gene));
      if (state.gene && genes.includes(state.gene)) $('color-mode').value = 'gene:' + state.gene;
      else { state.gene = ''; $('color-mode').value = 'state'; }
      const img = section.image || meta.image;
      $('image-control').hidden = !(img && img.path && Array.isArray(img.bounds) && img.bounds.length === 4);
      if (!$('image-control').hidden) { const image = new Image(); image.onload = () => { if (token === state.loadToken) { state.image = { image, bounds: img.bounds }; requestDraw(); } }; image.onerror = () => { if (token === state.loadToken) toast('The registered image is unavailable. Original spots remain visible.'); }; image.src = img.path; }
      refreshEdges(); fitView(); renderStateLegend();
    } catch (error) { if (token === state.loadToken) { setText('canvas-status', error.message); $('canvas-status').hidden = false; } }
  }

  function statCell(label, value, r, key) {
    const cell = node('div', 'stat-cell'); cell.append(node('label', null, label));
    const showValue = finite(value) && evaluable(r);
    const strong = node('strong', showValue ? '' : 'pending', showValue ? formatNumber(value, key) : missingLabel(r));
    strong.title = label + ' = ' + exact(value) + (!showValue && finite(value) ? ' · source placeholder' : ''); cell.append(strong); return cell;
  }
  function renderEvidence() {
    const unit = state.selected, r = currentResult(), stats = $('selected-stats'); empty(stats);
    renderDiagnostics(r);
    if (!unit) { setText('selected-lr', 'No candidate selected'); setText('selected-states', 'Adjust the filters to browse candidates.'); setText('selected-context', ''); setText('selected-status', ''); empty($('genetics-content')); empty($('donors-content')); return; }
    setText('selected-lr', candidateLabel(unit));
    setText('selected-states', unit.sender_state + ' → ' + unit.receiver_state);
    setText('selected-context', unit.context);
    stats.append(statCell(state.version.p_label || 'Candidate p', r && r.p, r, 'p'), statCell(correctionLabel(), r && r.q, r, 'q'));
    const score = r && finite(r.score) ? formatNumber(r.score) : missingLabel(r);
    setText('selected-status', 'Score ' + score + ' · ' + (r ? displayStatus(r.status, r) : 'No result in this version') + (!evaluable(r) && r && (finite(r.p) || finite(r.q)) ? '\nStored placeholders: p = ' + exact(r.p) + '; Adjusted value = ' + exact(r.q) : ''));
    $('selected-status').title = 'Hypothesis ' + unit.unit_id;
    if (r && r.evidence_scope) $('selected-status').append(document.createTextNode('\n' + r.evidence_scope));
    $('selected-status').style.whiteSpace = 'pre-line';
    renderGenetics(unit); renderDonors(unit);
  }
  function kv(items) {
    const list = node('dl', 'kv');
    for (const [label,value,title] of items) { list.append(node('dt', null, label)); const dd = node('dd', null, value == null || value === '' ? 'Not provided' : value); if (title) dd.title = title; list.append(dd); }
    return list;
  }
  function renderDiagnostics(result) {
    const container = $('selected-diagnostics'); empty(container);
    const keys = ['p_discrete', 'mcse_diagnostic', 'draws'];
    const present = result && keys.some(key => savedDiagnostic(result, key) !== null && savedDiagnostic(result, key) !== undefined);
    container.hidden = !present;
    if (!present) return;
    const raw = key => { const value = savedDiagnostic(result, key); return value == null || value === '' || /^(nan|null|none)$/i.test(String(value)) ? 'Not computed' : String(value); };
    const continuous = savedDiagnostic(result, 'p_continuous');
    const values = [['Model-tail p (continuous)', continuous == null ? exact(result.p) : String(continuous)], ['Empirical p (discrete)', raw('p_discrete')], ['MCSE diagnostic', raw('mcse_diagnostic')], ['Saved refits', raw('draws')]];
    container.append(node('summary', null, 'Tail diagnostics'), kv(values.map(([label,value]) => [label,value,'Source value: '+value])));
    if (!evaluable(result)) container.append(node('p', 'small-note', 'Source fields · unavailable hypothesis'));
  }
  function renderGenetics(unit) {
    const container = $('genetics-content'); empty(container);
    const links = report.genetic_links.filter(link => link.channel_id === unit.channel_id);
    if (!links.length) { container.append(node('div', 'empty-state', 'No genetic mapping records were provided for this channel.')); return; }
    const groups = new Map();
    for (const link of links) { const key = link.locus_id || 'Locus unspecified'; if (!groups.has(key)) groups.set(key, []); groups.get(key).push(link); }
    for (const [locus, values] of groups) {
      const group = node('section', 'genetic-group'); group.append(node('h3', null, locus));
      for (const link of values) {
        const card = node('div', 'evidence-card'), top = node('div', 'card-top');
        top.append(node('strong', null, link.rsid || link.variant_id || 'Variant unspecified'), node('span', 'tiny-badge', link.endpoint_role || 'Endpoint')); card.append(top);
        card.append(kv([['Variant / alleles', link.variant_id], ['rsID', link.rsid], ...(link.variant_role != null ? [['Variant role',link.variant_role]] : []), ['Gene', [link.gene_symbol, link.gene_id].filter(Boolean).join(' · ')], ['Mapping basis', link.mapping_basis], ...(link.gene_mapping_witness_variant_ids != null ? [['Mapping witnesses',Array.isArray(link.gene_mapping_witness_variant_ids) ? link.gene_mapping_witness_variant_ids.join('; ') : link.gene_mapping_witness_variant_ids]] : []), ['Historical PIP', finite(link.pip) ? formatNumber(link.pip, 'p') : 'Not provided', 'PIP = ' + exact(link.pip)], ['PIP status', link.pip_status], ['GWAS p', probability(link.gwas_p) ? formatNumber(link.gwas_p, 'p') : 'Not provided', 'GWAS p = ' + exact(link.gwas_p)], ['Credible set', link.credible_set], ['Genome build', link.build], ...(link.ld_resolution != null ? [['LD resolution', link.ld_resolution]] : [])]));
        if (link.mapping_note) { const details = node('details','source-note'); details.append(node('summary',null,'Mapping note'),node('p','small-note',link.mapping_note)); card.append(details); }
        group.append(card);
      }
      container.append(group);
    }
  }
  function renderDonors(unit) {
    const container = $('donors-content'); empty(container);
    if (sectionEvidence(metadata)) container.append(node('p','evidence-intro','Saved section evidence; no donor-level aggregation or replication test.'));
    const rows = report.donor_evidence.filter(row => row.unit_id === unit.unit_id && row.version === state.version.id);
    if (!rows.length) { container.append(node('div', 'empty-state', 'No evidence records for this version.')); return; }
    const maxScore = Math.max(...rows.map(row => finite(row.score) ? Math.abs(row.score) : 0), 1e-12);
    for (const row of rows) {
      const card = node('div', 'evidence-card'), top = node('div', 'card-top'); top.append(node('strong', null, row.donor_id), node('span', 'tiny-badge', displayStatus(row.status,row))); card.append(top);
      if (finite(row.score)) { const bar = node('div', 'donor-spark'); bar.title = 'Score: ' + exact(row.score); const fill = node('i'); fill.style.width = Math.abs(row.score) / maxScore * 100 + '%'; if (row.score < 0) fill.style.background = '#bda285'; bar.append(fill); card.append(bar); }
      const placeholder = unavailableStatus(row.status) || row.evaluable === false;
      card.append(kv([[sectionEvidence(metadata) ? 'Section score' : 'Donor score', finite(row.score) ? formatNumber(row.score) : 'Not computed', 'score = ' + exact(row.score)], [sectionEvidence(metadata) ? 'Section p' : 'Donor p', probability(row.p) ? formatNumber(row.p, 'p') + (placeholder ? ' · source placeholder' : '') : placeholder ? 'Not evaluable' : 'Not computed', 'Source p = ' + exact(row.p)], ['Sections', Array.isArray(row.sections) ? row.sections.length : row.sections], ['Spatial edges', row.edges]]));
      if (row.evidence_scope && metadata.analysis_scope !== 'single_embryo_section') { const details = node('details','source-note'); details.append(node('summary',null,'Evidence scope'),node('p','small-note',row.evidence_scope)); card.append(details); }
      if (row.loo) { card.append(node('p', 'small-note', 'Saved leave-one-donor-out (LOO) results'), kv([['LOO score', finite(row.loo.score) ? formatNumber(row.loo.score) : 'Not computed'], ['LOO p', probability(row.loo.p) ? formatNumber(row.loo.p, 'p') : 'Not computed', 'LOO p = ' + exact(row.loo.p)]])); }
      container.append(card);
    }
    container.append(node('p', 'small-note', 'Bars: relative score magnitude.'));
  }

  function refreshEdges() {
    state.edges = state.section && state.sectionMeta && state.selected && state.sectionMeta.context === state.selected.context ? orientedEdges(state.section, state.selected) : [];
    requestDraw(); updatePlotSummary();
  }
  function updatePlotSummary() {
    if (!state.section) return;
    const r = currentResult(), adj = adjustment(state.version);
    const sig = state.highlight && adj.eligible && evaluable(r) && probability(r.q) && r.q < .05;
    setText('plot-summary', state.edges.length.toLocaleString() + ' matching spatial edges · ' + (evaluable(r) ? 'Candidate p = ' + formatNumber(r.p, 'p') : missingLabel(r) + ' (gray)') + (sig ? '\n' + adj.label + ' < 0.05 · dark outline' : '') + (state.gene ? '\nExpression: ' + state.gene + ' · ' + (state.section.expression_label || state.sectionMeta.expression_label || 'Source values') : ''));
    $('plot-summary').style.whiteSpace = 'pre-line';
    $('plot-summary').title = 'Candidate p = ' + exact(r && r.p) + '; ' + correctionLabel() + ' = ' + exact(r && r.q) + '\nColor domain: 0–6.';
  }
  function expressionDomain() {
    if (!state.gene || !state.section) return null;
    const values = state.section.expression[state.gene] || [];
    let min = Infinity, max = -Infinity;
    for (const value of values) if (finite(value)) { min = Math.min(min, value); max = Math.max(max, value); }
    return min === Infinity ? null : { min, max };
  }
  function expressionColor(value, domain) {
    if (!finite(value) || !domain) return '#d8ded7';
    const t = domain.max === domain.min ? .5 : (value - domain.min) / (domain.max - domain.min);
    return 'rgb(' + [229,237,204].map((v,i) => Math.round(v + t * ([29,115,106][i] - v))).join(',') + ')';
  }
  function renderStateLegend() {
    const legend = $('state-legend'); empty(legend);
    if (!state.section) return;
    if (state.gene) { const d = expressionDomain(); legend.append(node('span', null, d ? state.gene + ': ' + formatNumber(d.min) + ' (light) → ' + formatNumber(d.max) + ' (dark); missing values are gray' : state.gene + ': no valid expression values')); return; }
    for (const name of unique(state.section.spots.map(s => s.state))) { const item = node('span'), dot = node('i'); dot.style.background = stateColor(name); item.append(dot, document.createTextNode(name)); legend.append(item); }
  }
  function resizeCanvas() {
    const wrap = $('canvas-wrap');
    if (wrap.clientWidth < 90 || wrap.clientHeight < 90) return;
    width = wrap.clientWidth; height = wrap.clientHeight;
    const dpr = Math.min(root.devicePixelRatio || 1, 3);
    if (canvas.width !== Math.round(width * dpr)) canvas.width = Math.round(width * dpr);
    if (canvas.height !== Math.round(height * dpr)) canvas.height = Math.round(height * dpr);
    if (state.view && state.bounds) state.view = resizedView(state.view, fitScaleForBounds(state.bounds, width, height));
    else if (state.section) fitView();
    requestDraw();
  }
  function orientedY(y) { return state.sectionMeta && state.sectionMeta.coordinate_system && state.sectionMeta.coordinate_system.y_direction === 'up' ? -y : y; }
  function fitView() {
    if (!state.section || !state.section.spots.length) return;
    let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;
    for (const spot of state.section.spots) { xmin=Math.min(xmin,spot.x); xmax=Math.max(xmax,spot.x); ymin=Math.min(ymin,orientedY(spot.y)); ymax=Math.max(ymax,orientedY(spot.y)); }
    state.bounds = { xmin, xmax, ymin, ymax };
    const scale = fitScaleForBounds(state.bounds, width, height);
    state.view = { cx:(xmin+xmax)/2,cy:(ymin+ymax)/2,scale,fitScale:scale }; requestDraw();
  }
  function project(spot) { return { x:width/2+(spot.x-state.view.cx)*state.view.scale,y:height/2+(orientedY(spot.y)-state.view.cy)*state.view.scale }; }
  function unproject(x,y) { return { x:state.view.cx+(x-width/2)/state.view.scale, y:state.view.cy+(y-height/2)/state.view.scale }; }
  function requestDraw() { if (drawPending) return; drawPending = true; requestAnimationFrame(() => { drawPending = false; draw(); }); }
  function draw() {
    const dpr = canvas.width / width;
    ctx.setTransform(dpr,0,0,dpr,0,0); ctx.clearRect(0,0,width,height); ctx.fillStyle='#fbfcf9'; ctx.fillRect(0,0,width,height);
    if (!state.section || !state.view) return;
    const section=state.section, view=state.view, r=currentResult(), adj=adjustment(state.version);
    const highlighted=state.highlight && adj.eligible && evaluable(r) && probability(r.q) && r.q<.05;
    const score=evaluable(r) ? Math.min(6,r.p===0?6:-Math.log10(r.p)) : null;
    const lineWidth=score===null?1:1+score/6*2.5;
    ctx.save(); ctx.beginPath(); ctx.rect(35,14,width-49,height-46); ctx.clip();
    if (state.showImage && state.image) {
      const bounds=state.image.bounds, a=project({x:bounds[0],y:bounds[1]}),b=project({x:bounds[2],y:bounds[3]});
      ctx.globalAlpha=.52;ctx.drawImage(state.image.image,Math.min(a.x,b.x),Math.min(a.y,b.y),Math.abs(b.x-a.x),Math.abs(b.y-a.y));ctx.globalAlpha=1;
    }
    const positions=section.spots.map(project), radius=Math.max(1.45,Math.min(4.6,1.7*Math.sqrt(view.scale/view.fitScale))),domain=expressionDomain();
    for (let i=0;i<section.spots.length;i++) {
      const point=positions[i]; if(point.x<30||point.x>width||point.y<0||point.y>height)continue;
      const spot=section.spots[i], relevant=state.selected && (spot.state===state.selected.sender_state||spot.state===state.selected.receiver_state);
      ctx.globalAlpha = state.gene ? .88 : relevant ? .8 : .44;
      ctx.fillStyle=state.gene?expressionColor(section.expression[state.gene][i],domain):stateColor(spot.state);
      ctx.beginPath();ctx.arc(point.x,point.y,radius,0,Math.PI*2);ctx.fill();
    }
    ctx.globalAlpha=1;
    if (state.showEdges && state.edges.length) {
      const edgeColor=pColor(evaluable(r) ? r.p : null);
      function drawEdges(outline) {
        ctx.strokeStyle=outline?'#28483b':edgeColor;ctx.fillStyle=outline?'#28483b':edgeColor;ctx.lineWidth=outline?lineWidth+2:lineWidth;ctx.globalAlpha=outline ? .9 : .83;ctx.lineCap='round';
        ctx.beginPath();
        for (const edge of state.edges) { const a=positions[edge[0]],b=positions[edge[1]];if(Math.max(a.x,b.x)<30||Math.min(a.x,b.x)>width||Math.max(a.y,b.y)<0||Math.min(a.y,b.y)>height)continue;ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y); }
        ctx.stroke();
      }
      if (highlighted)drawEdges(true);drawEdges(false);
      ctx.globalAlpha=.9;
      for (const edge of state.edges) {
        const a=positions[edge[0]],b=positions[edge[1]],dx=b.x-a.x,dy=b.y-a.y,len=Math.hypot(dx,dy);
        if(len<10||Math.max(a.x,b.x)<30||Math.min(a.x,b.x)>width||Math.max(a.y,b.y)<0||Math.min(a.y,b.y)>height)continue;
        const ux=dx/len,uy=dy/len,arrow=Math.min(5,lineWidth+2),x=a.x+dx*.69,y=a.y+dy*.69;
        ctx.beginPath();ctx.moveTo(x+ux*arrow,y+uy*arrow);ctx.lineTo(x-ux*arrow-uy*arrow*.65,y-uy*arrow+ux*arrow*.65);ctx.lineTo(x-ux*arrow+uy*arrow*.65,y-uy*arrow-ux*arrow*.65);ctx.closePath();ctx.fill();
      }
      ctx.globalAlpha=1;
    }
    if (state.dragging && state.mode==='box') { const d=state.dragging;ctx.fillStyle='#17695f16';ctx.strokeStyle='#17695f';ctx.lineWidth=1;ctx.setLineDash([5,3]);ctx.fillRect(d.startX,d.startY,d.x-d.startX,d.y-d.startY);ctx.strokeRect(d.startX,d.startY,d.x-d.startX,d.y-d.startY);ctx.setLineDash([]); }
    if (state.hovered!==null) { const p=positions[state.hovered];if(p){ctx.strokeStyle='#244c3e';ctx.lineWidth=1.5;ctx.beginPath();ctx.arc(p.x,p.y,radius+3,0,Math.PI*2);ctx.stroke();} }
    ctx.restore(); drawAxes();
  }
  function drawAxes() {
    ctx.strokeStyle='#dfe6dd';ctx.fillStyle='#8b9a8a';ctx.lineWidth=.7;ctx.font='9px Arial,sans-serif';ctx.textAlign='center';
    ctx.beginPath();ctx.moveTo(35,14);ctx.lineTo(35,height-32);ctx.lineTo(width-14,height-32);ctx.stroke();
    for(let i=0;i<=4;i++) {
      const x=35+(width-49)*i/4,y=14+(height-46)*i/4,worldX=unproject(x,0).x,worldY=orientedY(unproject(0,y).y);
      ctx.fillText(axisNumber(worldX),x,height-17);ctx.save();ctx.textAlign='right';ctx.fillText(axisNumber(worldY),30,y+3);ctx.restore();
    }
    ctx.textAlign='right';ctx.fillText('x',width-6,height-18);ctx.textAlign='left';ctx.fillText('y',9,11);
  }
  function axisNumber(value) { if(Math.abs(value)>=1e5 || Math.abs(value)>0 && Math.abs(value)<.01)return value.toExponential(1);return Number(value.toPrecision(3)).toString(); }
  function pointer(e) { const rect=canvas.getBoundingClientRect();return {x:e.clientX-rect.left,y:e.clientY-rect.top}; }
  function startPointer(e) {
    if (!state.view || e.button!==0) return;
    e.preventDefault();canvas.focus({preventScroll:true});canvas.setPointerCapture(e.pointerId);
    const p=pointer(e);state.dragging={startX:p.x,startY:p.y,x:p.x,y:p.y,cx:state.view.cx,cy:state.view.cy};state.hovered=null;$('canvas-tooltip').hidden=true;canvas.classList.add('is-dragging');requestDraw();
  }
  function movePointer(e) {
    if(!state.view)return;
    const p=pointer(e);
    if(state.dragging) { const d=state.dragging;d.x=p.x;d.y=p.y;if(state.mode==='pan'){state.view.cx=d.cx-(p.x-d.startX)/state.view.scale;state.view.cy=d.cy-(p.y-d.startY)/state.view.scale;}requestDraw();return; }
    let nearest=null,best=64;
    state.section.spots.forEach((spot,i)=>{const pos=project(spot),distance=(pos.x-p.x)**2+(pos.y-p.y)**2;if(distance<best){nearest=i;best=distance;}});
    if(nearest!==state.hovered){state.hovered=nearest;requestDraw();}
    const tip=$('canvas-tooltip');
    if(nearest===null){tip.hidden=true;return;}
    const spot=state.section.spots[nearest];
    tip.textContent=[currentPointKind().singular+' '+spot.id,spot.state,'x = '+exact(spot.x)+' · y = '+exact(spot.y),state.gene?state.gene+' = '+exact(state.section.expression[state.gene][nearest]):'',currentPointKind().description].filter(Boolean).join('\n');tip.style.whiteSpace='pre-line';tip.hidden=false;tip.style.left=Math.max(5,Math.min(width-250,p.x+14))+'px';tip.style.top=Math.max(5,Math.min(height-135,p.y+12))+'px';
  }
  function endPointer(e) {
    if(!state.dragging)return;
    const d=state.dragging;
    if(state.mode==='box'&&Math.abs(d.x-d.startX)>10&&Math.abs(d.y-d.startY)>10){const center=unproject((d.x+d.startX)/2,(d.y+d.startY)/2),factor=Math.min((width-75)/Math.abs(d.x-d.startX),(height-65)/Math.abs(d.y-d.startY));state.view.cx=center.x;state.view.cy=center.y;state.view.scale=Math.min(state.view.fitScale*100,state.view.scale*factor);toast('Region selected.');}
    state.dragging=null;canvas.classList.remove('is-dragging');if(canvas.hasPointerCapture(e.pointerId))canvas.releasePointerCapture(e.pointerId);requestDraw();
  }
  function zoom(e) {
    if(!state.view)return;e.preventDefault();const p=pointer(e),before=unproject(p.x,p.y),factor=Math.exp(-Math.max(-180,Math.min(180,e.deltaY))*.0018);
    state.view.scale=Math.max(state.view.fitScale*.15,Math.min(state.view.fitScale*100,state.view.scale*factor));const after=unproject(p.x,p.y);state.view.cx+=before.x-after.x;state.view.cy+=before.y-after.y;state.hovered=null;$('canvas-tooltip').hidden=true;requestDraw();
  }
  function saveBlob(blob,filename) { const url=URL.createObjectURL(blob),a=node('a');a.href=url;a.download=filename;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),20000); }
  function safeName(value) { return String(value || 'report').replace(/[^a-zA-Z0-9._-]/g,'_').slice(0,110); }
  function exportPNG() {
    if(!state.section){toast('No section is available for export.');return;}
    state.hovered=null;$('canvas-tooltip').hidden=true;canvas.width=canvas.width;draw();
    const output=document.createElement('canvas'),scale=2,w=Math.max(820,Math.round(width)),plotHeight=Math.round(height*w/width),header=106,footer=130;
    output.width=w*scale;output.height=(plotHeight+header+footer)*scale;
    const g=output.getContext('2d');g.scale(scale,scale);g.fillStyle='#fff';g.fillRect(0,0,w,output.height/scale);g.fillStyle='#17695f';g.font='bold 22px Arial,sans-serif';g.fillText('LUCENT  /  '+(state.selected?candidateLabel(state.selected):'Spatial section'),24,34,w-48);
    g.fillStyle='#6c7c72';g.font='12px Arial,sans-serif';g.fillText((state.sectionMeta.label||state.sectionMeta.id)+' · '+state.sectionMeta.context+' · '+replicateLabel(metadata)+' '+state.sectionMeta.donor_id+' · '+currentPointKind().description,24,58,w-48);g.fillText('Version: '+state.version.id+' | '+(state.version.label||'')+' | '+(state.version.adjustment_method||'adjustment unavailable'),24,79,w-48);
    g.drawImage(canvas,0,header,w,plotHeight);
    const r=currentResult(),y=header+plotHeight+23;
    g.fillStyle='#3c5146';g.fillText((evaluable(r)?'Candidate-level ':r&&finite(r.p)?'Source placeholder · ':missingLabel(r)+' · ')+(state.version.p_label||'p')+' = '+exact(r&&r.p)+'  |  '+correctionLabel()+' = '+exact(r&&r.q)+'  |  '+state.edges.length+' matching spatial edges',24,y,w-48);
    const gradient=g.createLinearGradient(24,0,220,0);ramp.forEach((c,i)=>gradient.addColorStop(i/(ramp.length-1),'rgb('+c.join(',')+')'));g.fillStyle=gradient;g.fillRect(24,y+14,196,8);g.fillStyle='#6c7c72';g.font='10px Arial,sans-serif';g.fillText('−log₁₀(p):  0                      3                       ≥6',24,y+37);
    g.fillText(state.gene?'Expression: '+state.gene+' · '+(state.section.expression_label||state.sectionMeta.expression_label||'source values'):'States: '+unique(state.section.spots.map(s=>s.state)).join(' · '),245,y+22,w-269);
    g.fillText('Color and width: candidate-level p. Arrows: ligand → receptor. Equal-aspect coordinates.',24,y+61,w-48);
    g.fillText('Highlight: '+(state.highlight?adjustment(state.version).highlight+' < 0.05':'Off')+' · Edges: '+(state.showEdges?'Shown':'Hidden')+' · Unit: '+((state.sectionMeta.coordinate_system||{}).unit||'source coordinate'),24,y+80,w-48);
    try{output.toBlob(blob=>{if(blob){saveBlob(blob,'LUCENT_'+safeName(state.version.id)+'_'+safeName(state.sectionMeta.id)+'.png');toast('PNG exported.');}},'image/png');}catch(error){toast('Image export blocked. Turn off the tissue image and retry.');}
  }
  function bindEvents() {
    $('version').addEventListener('change',()=>{state.version=versions.find(v=>v.id===$('version').value);state.filters.evidence=report.units.some(u=>evaluable(resultFor(u,state.version.id)))?'evaluable':'all';state.page=0;updateVersion();updateCandidates();});
    for(const id of ['context','sender','receiver'])$(id).addEventListener('change',()=>{state.filters[id]=$(id).value;state.page=0;if(id==='context')updateSections();updateCandidates();});
    $('search').addEventListener('input',()=>{state.filters.search=$('search').value;state.page=0;updateCandidates();});
    document.querySelectorAll('input[name=evidence]').forEach(el=>el.addEventListener('change',()=>{state.filters.evidence=el.value;state.page=0;updateCandidates();}));
    $('donor').addEventListener('change',()=>updateSections());
    $('section').addEventListener('change',()=>loadSection(report.sections.find(s=>s.id===$('section').value)||null));
    $('clear-filters').addEventListener('click',()=>{state.filters={context:'',sender:'',receiver:'',search:'',evidence:report.units.some(u=>evaluable(resultFor(u,state.version.id)))?'evaluable':'all'};for(const id of ['context','sender','receiver','search'])$(id).value='';$('donor').value='';state.page=0;updateVersion();updateSections(metadata.default_section_id);updateCandidates(false);});
    for(const mode of ['pan','box'])$('mode-'+mode).addEventListener('click',()=>{state.mode=mode;for(const name of ['pan','box']){$('mode-'+name).classList.toggle('active',name===mode);$('mode-'+name).setAttribute('aria-pressed',String(name===mode));}canvas.classList.toggle('box-mode',mode==='box');setText('interaction-hint',mode==='box'?'Drag a rectangle to zoom':'Scroll to zoom · drag to pan');});
    $('reset-view').addEventListener('click',fitView);
    $('color-mode').addEventListener('change',()=>{state.gene=$('color-mode').value.startsWith('gene:')?$('color-mode').value.slice(5):'';renderStateLegend();updatePlotSummary();requestDraw();});
    $('show-edges').addEventListener('change',()=>{state.showEdges=$('show-edges').checked;requestDraw();});
    $('highlight-adjusted').addEventListener('change',()=>{state.highlight=$('highlight-adjusted').checked;updatePlotSummary();requestDraw();});
    $('show-image').addEventListener('change',()=>{state.showImage=$('show-image').checked;requestDraw();});
    for(const name of ['genetics','donors'])$('tab-'+name).addEventListener('click',()=>{for(const tab of ['genetics','donors']){$('tab-'+tab).classList.toggle('active',tab===name);$('tab-'+tab).setAttribute('aria-selected',String(tab===name));$(tab+'-content').hidden=tab!==name;}});
    $('page-prev').addEventListener('click',()=>{state.page--;renderTable();document.querySelector('.table-scroll').scrollTop=0;});
    $('page-next').addEventListener('click',()=>{state.page++;renderTable();document.querySelector('.table-scroll').scrollTop=0;});
    document.querySelectorAll('[data-sort]').forEach(button=>button.addEventListener('click',()=>{const key=button.dataset.sort;state.ascending=state.sort===key?!state.ascending:key!=='score';state.sort=key;state.page=0;updateCandidates();}));
    $('export-tsv').addEventListener('click',()=>{saveBlob(new Blob(['\ufeff'+unitsTSV(state.filtered,state.version)],{type:'text/tab-separated-values;charset=utf-8'}),'LUCENT_'+safeName(state.version.id)+'_filtered_candidates.tsv');toast(state.filtered.length+' candidates exported.');});
    $('export-png').addEventListener('click',exportPNG);
    canvas.addEventListener('pointerdown',startPointer);canvas.addEventListener('pointermove',movePointer);canvas.addEventListener('pointerup',endPointer);canvas.addEventListener('pointercancel',()=>{state.dragging=null;canvas.classList.remove('is-dragging');requestDraw();});canvas.addEventListener('pointerleave',()=>{if(!state.dragging){state.hovered=null;$('canvas-tooltip').hidden=true;requestDraw();}});canvas.addEventListener('wheel',zoom,{passive:false});
    canvas.addEventListener('keydown',e=>{if(!state.view)return;const step=40/state.view.scale;if(['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','+','-','0','Escape'].includes(e.key))e.preventDefault();if(e.key==='ArrowLeft')state.view.cx-=step;if(e.key==='ArrowRight')state.view.cx+=step;if(e.key==='ArrowUp')state.view.cy-=step;if(e.key==='ArrowDown')state.view.cy+=step;if(e.key==='+')state.view.scale=Math.min(state.view.fitScale*100,state.view.scale*1.2);if(e.key==='-')state.view.scale=Math.max(state.view.fitScale*.15,state.view.scale/1.2);if(e.key==='0')fitView();if(e.key==='Escape'){state.dragging=null;$('canvas-tooltip').hidden=true;}requestDraw();});
  }
  // Read-only handle for local QA; changing filters still goes through the interface.
  root.LUCENT_VIEW = { getState: () => ({ version:state.version.id,selected:state.selected&&state.selected.unit_id,section:state.sectionMeta&&state.sectionMeta.id,filtered:state.filtered.length,edges:state.edges.length,view:state.view?{...state.view}:null }), redraw:requestDraw };
  initialize();
})(typeof window !== 'undefined' ? window : globalThis);
