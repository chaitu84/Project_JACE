"use strict";
(() => {
  const $ = id => document.getElementById(id);
  const escape = value => String(value).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const fmt = value => Number(value).toLocaleString('en-US');
  let preview = [];
  let targetLength = 0;
  let selectedMatch = null;
  const k = $('rna-k');
  k.addEventListener('input', () => { $('rna-k-value').textContent = `${k.value} nt`; });

  async function loadData() {
    try {
      const [statusResponse, targetResponse] = await Promise.all([fetch('/api/status'), fetch('/api/target')]);
      if (!statusResponse.ok) throw new Error('Local analysis server is unavailable');
      const status = await statusResponse.json();
      $('rna-availability').textContent = status.target_available && status.human_available
        ? `Ready · ${status.human_source}` : 'Input data missing';
      if (targetResponse.ok) {
        const target = await targetResponse.json();
        targetLength = target.length_nt;
        $('rna-target-summary').textContent = `${fmt(target.length_nt)} nucleotides · ${target.id} · exact sequence from local FASTA`;
        $('rna-target-end').textContent = `3′ · ${fmt(target.length_nt)} nt`;
        $('rna-coverage-end').textContent = `3′ · ${fmt(target.length_nt)} nt`;
        $('rna-target').innerHTML = `<p>${escape(target.id)} · ${fmt(target.length_nt)} nt · spliced viral target</p><p class="rna-target-sequence">${escape(target.sequence)}</p>`;
      } else $('rna-target').textContent = 'Viral target FASTA not found on server.';
    } catch (error) {
      $('rna-availability').textContent = 'Start the Flask server to enable live analysis';
      $('rna-target').textContent = 'Target sequence is available when Flask is running.';
    }
  }
  function bars(container, entries) {
    const max = Math.max(1, ...entries.map(entry => Number(entry[1])));
    container.innerHTML = entries.length ? entries.map(([label, count]) =>
      `<div class="rna-bar-row"><span title="${escape(label)}">${escape(label)}</span><div class="rna-bar-track"><div class="rna-bar-fill" style="width:${100*count/max}%"></div></div><strong>${fmt(count)}</strong></div>`).join('') : '<p>No matches.</p>';
  }
  function baseLine(sequence) {
    return [...sequence].map(base => `<span class="rna-base rna-base-${escape(base)}">${escape(base)}</span>`).join('');
  }
  function showAlignment(index) {
    const match = preview[index];
    if (!match) return;
    const viral = match['Viral Sequence'];
    const human = match['Human Matched Sequence'];
    $('rna-alignment').innerHTML = `<h3>Antiparallel alignment · ${escape(match['Human Transcript'])}</h3>
      <p>Human ${match['Human Start']}–${match['Human End']} · Viral ${match['Viral Start']}–${match['Viral End']} · zero-based, end-exclusive</p>
      <pre>Viral 5′ ${baseLine(viral)} 3′\nHuman 3′ ${baseLine([...human].reverse().join(''))} 5′</pre>
      <p class="rna-footnote">Exact reverse complement; mismatches 0, gaps 0. This is a sequence screen, not an energy prediction.</p>`;
    selectedMatch = match;
    highlightSelectedMatch();
    $('rna-alignment').scrollIntoView({behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'instant' : 'smooth', block: 'center'});
  }
  function renderTable() {
    const query = $('rna-filter').value.trim().toLowerCase();
    const matches = preview.map((match, index) => ({match, index})).filter(({match}) =>
      match['Human Transcript'].toLowerCase().includes(query));
    $('rna-table-body').innerHTML = matches.map(({match, index}) => `<tr>
      <td>${escape(match['Human Transcript'])}</td><td>${match['Human Start']}–${match['Human End']}</td>
      <td>${match['Viral Start']}–${match['Viral End']}</td><td>${match['Match Length']} nt</td>
      <td><button type="button" data-rna-index="${index}">View alignment</button></td></tr>`).join('') || '<tr><td colspan="5">No matching preview rows.</td></tr>';
  }
  $('rna-filter').addEventListener('input', renderTable);
  $('rna-table-body').addEventListener('click', event => {
    const button = event.target.closest('[data-rna-index]');
    if (button) showAlignment(Number(button.dataset.rnaIndex));
  });
  $('rna-form').addEventListener('submit', async event => {
    event.preventDefault();
    const submit = $('rna-submit');
    submit.disabled = true;
    $('rna-status').textContent = 'Scanning human transcripts. Keep this page open…';
    try {
      const response = await fetch('/api/analyze', {method:'POST', headers:{'Content-Type':'application/json'},
        body:JSON.stringify({match_length:Number(k.value), transcript_limit:Number($('rna-limit').value)})});
      const data = await response.json();
      if (!response.ok) throw new Error(data.error || 'Analysis failed');
      preview = data.preview;
      selectedMatch = null;
      $('rna-status').textContent = `Completed · ${fmt(data.analyzed)} transcript records scanned · ${fmt(data.total_matches)} exact matches.`;
      $('rna-results').hidden = false;
      $('rna-metrics').innerHTML = [
        ['Exact matches', data.total_matches], ['Unique matching transcripts', data.unique_transcripts],
        ['Records scanned', data.analyzed], ['Match length (nt)', data.match_length]
      ].map(([label,value]) => `<div class="rna-metric"><strong>${fmt(value)}</strong><span>${label}</span></div>`).join('');
      $('rna-download').href = `/api/download/${encodeURIComponent(data.run_id)}`;
      $('rna-filter').value = '';
      renderTable();
      $('rna-alignment').innerHTML = '<h3>Selected transcript alignment</h3><p>Select “View alignment” in the results table below to inspect a real match here.</p>';
      bars($('rna-transcript-chart'), data.transcript_counts);
      renderCoverage(data);
      const counts = new Map();
      const binSize = Math.max(1, Math.ceil(data.viral_length / 20));
      for (const [position, count] of data.position_counts) {
        const start = Math.floor(position / binSize) * binSize;
        counts.set(start, (counts.get(start) || 0) + count);
      }
      bars($('rna-position-chart'), [...counts].sort((a,b) => a[0]-b[0]).map(([start,count]) => [`${start}–${Math.min(data.viral_length-1,start+binSize-1)}`,count]));
    } catch (error) { $('rna-status').textContent = `Analysis error: ${error.message}`; }
    finally { submit.disabled = false; }
  });
  function renderCoverage(data) {
    const length = data.viral_length;
    targetLength = length;
    const cells = 48;
    const totals = Array(cells).fill(0);
    for (const [start, count] of data.position_counts) {
      const first = Math.min(cells - 1, Math.floor(start / length * cells));
      const last = Math.min(cells - 1, Math.floor((start + data.match_length - 1) / length * cells));
      for (let i = first; i <= last; i++) totals[i] += count;
    }
    const peak = Math.max(1, ...totals);
    $('rna-coverage').innerHTML = totals.map((value, i) => {
      const from = Math.floor(i * length / cells);
      const to = Math.min(length - 1, Math.ceil((i+1)*length/cells)-1);
      const opacity = value ? (0.18 + 0.82 * Math.sqrt(value / peak)).toFixed(3) : '0.06';
      return `<div class="rna-coverage-cell" data-from="${from}" data-to="${to}" tabindex="0" style="--coverage-opacity:${opacity}" title="Target positions ${from}–${to}: ${fmt(value)} overlapping matches" aria-label="Positions ${from} to ${to}, ${value} overlapping matches"></div>`;
    }).join('');
    $('rna-coverage-summary').textContent = `${fmt(data.total_matches)} matches · ${fmt(length)}-nt target`;
    $('rna-coverage-end').textContent = `3′ · ${fmt(length)} nt`;
    highlightSelectedMatch();
  }
  function highlightSelectedMatch() {
    const cells = $('rna-coverage').querySelectorAll('.rna-coverage-cell');
    for (const cell of cells) {
      const from = Number(cell.dataset.from);
      const to = Number(cell.dataset.to);
      const active = selectedMatch && from < Number(selectedMatch['Viral End']) && to >= Number(selectedMatch['Viral Start']);
      cell.classList.toggle('rna-coverage-selected', Boolean(active));
      if (active) {
        cell.setAttribute('aria-current', 'location');
      } else {
        cell.removeAttribute('aria-current');
      }
    }
    const summary = $('rna-coverage-summary');
    if (selectedMatch) {
      summary.textContent = `Selected viral span: ${fmt(selectedMatch['Viral Start'])}–${fmt(selectedMatch['Viral End'])} · ${fmt(targetLength)}-nt target`;
    }
  }
  function drawHelix() {
    const svg = $('rna-helix-svg');
    if (!svg) return;
    const pieces = [];
    for (let i = 0; i <= 30; i++) {
      const x = 22 + i * 11.8;
      const phase = i * Math.PI / 7;
      const y1 = 105 + 60 * Math.sin(phase);
      const y2 = 105 - 60 * Math.sin(phase);
      const opacity = 0.3 + 0.7 * Math.abs(Math.cos(phase));
      pieces.push(`<line x1="${x}" y1="${y1}" x2="${x}" y2="${y2}" stroke="#9cbeb4" stroke-width="2.4" opacity="${opacity}"/>`);
      pieces.push(`<circle cx="${x}" cy="${y1}" r="4" fill="#0e6964"/><circle cx="${x}" cy="${y2}" r="4" fill="#bd684e"/>`);
    }
    svg.innerHTML = `<g class="rna-helix-motion">${pieces.join('')}</g>`;
  }
  drawHelix();
  loadData();
})();
