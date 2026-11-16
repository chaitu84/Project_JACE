/* Charts and interactions use only the verified, saved control energies. */
"use strict";

const COLORS = { ink: "#172d32", muted: "#647674", teal: "#216f68", bar: "#5b9387", grid: "#e1e6de", original: "#b54e38", paper: "#fcfcf8" };
const number = new Intl.NumberFormat("en-US");
const energy = value => Number(value).toFixed(2);
const fraction = value => Number(value).toFixed(4);
const escapeHTML = value => String(value).replace(/[&<>"']/g, character => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[character]));

function histogram(values, edges) {
  const bins = edges.slice(0, -1).map((start, index) => ({ start, end: edges[index + 1], count: 0, entries: [] }));
  values.forEach((value, index) => {
    if (!Number.isFinite(value)) throw new Error("A control energy is missing or nonfinite.");
    const bin = bins.find((item, binIndex) => value >= item.start && (value < item.end || (binIndex === bins.length - 1 && value <= item.end)));
    if (!bin) throw new Error("A control energy lies outside the shared chart bins.");
    bin.count += 1;
    bin.entries.push({ shuffle: index + 1, energy: value });
  });
  return bins;
}

function chartSVG(candidate, edges, axisMax, idPrefix) {
  const width = 420, height = 260;
  const left = 44, right = 14, top = 38, bottom = 48;
  const plotWidth = width - left - right, plotHeight = height - top - bottom;
  const low = edges[0], high = edges[edges.length - 1];
  const x = value => left + (value - low) / (high - low) * plotWidth;
  const y = value => top + plotHeight - value / axisMax * plotHeight;
  const bins = histogram(candidate.energies, edges);
  const id = `${idPrefix}-${candidate.id.replace(/[^a-zA-Z0-9]/g, "-")}`;
  const originalX = x(candidate.original_energy);
  let ticks = "";
  const countStep = Math.max(1, Math.ceil(axisMax / 5));
  for (let count = 0; count <= axisMax; count += countStep) {
    ticks += `<line x1="${left}" y1="${y(count)}" x2="${width - right}" y2="${y(count)}" stroke="${COLORS.grid}"/>`;
    ticks += `<text x="${left - 9}" y="${y(count) + 4}" text-anchor="end" fill="${COLORS.muted}" font-size="14">${count}</text>`;
  }
  const tickSpacing = Math.max(2, Math.ceil((high - low) / 7 / 2) * 2);
  for (let tick = low; tick <= high; tick += tickSpacing) {
    ticks += `<text x="${x(tick)}" y="${height - bottom + 20}" text-anchor="middle" fill="${COLORS.muted}" font-size="14">${tick}</text>`;
  }
  const bars = bins.map((bin, index) => {
    const barWidth = x(bin.end) - x(bin.start) - 2;
    return `<rect class="hist-bar" data-bin="${index}" x="${x(bin.start) + 1}" y="${y(bin.count)}" width="${barWidth}" height="${plotHeight * bin.count / axisMax}" rx="1.5" fill="${COLORS.bar}"><title>${bin.start} to ${bin.end} kcal/mol, ${bin.count} controls</title></rect>`;
  }).join("");
  return `<svg class="histogram" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" tabindex="0" aria-labelledby="${id}-title ${id}-desc" data-candidate="${escapeHTML(candidate.id)}">
    <title id="${id}-title">${escapeHTML(candidate.id)} shuffled interaction energies</title>
    <desc id="${id}-desc">100 actual controls. Original energy ${energy(candidate.original_energy)} kcal/mol. ${candidate.controls_le_original} controls have energy at or below the original. All five plots share bins and axes. Use left and right arrow keys to inspect bins.</desc>
    <g font-family="Arial, sans-serif">${ticks}${bars}
    <line x1="${left}" y1="${y(0)}" x2="${width - right}" y2="${y(0)}" stroke="${COLORS.muted}" stroke-width="0.8"/>
    <line x1="${originalX}" y1="${top - 3}" x2="${originalX}" y2="${y(0)}" stroke="${COLORS.original}" stroke-width="2.5"/>
    <circle cx="${originalX}" cy="${top - 4}" r="3.5" fill="${COLORS.original}"/>
    <text x="${originalX}" y="${top - 16}" text-anchor="middle" fill="${COLORS.original}" font-size="14" font-weight="600">Original ${energy(candidate.original_energy)}</text>
    <text x="${left + plotWidth / 2}" y="${height - 6}" text-anchor="middle" fill="${COLORS.muted}" font-size="13">Interaction energy (kcal/mol)</text>
    <text x="13" y="${top + plotHeight / 2}" transform="rotate(-90 13 ${top + plotHeight / 2})" text-anchor="middle" fill="${COLORS.muted}" font-size="13">Control count</text>
    </g></svg>`;
}

function geneLabel(candidate) {
  return candidate.gene_name && !candidate.gene_name.startsWith("ENSG") ? candidate.gene_name : candidate.gene_id;
}

function downloadFile(filename, contents, type) {
  const url = URL.createObjectURL(new Blob([contents], { type }));
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function attachChartTooltip(wrapper, candidate, edges) {
  const svg = wrapper.querySelector("svg");
  const bins = histogram(candidate.energies, edges);
  const tooltip = document.createElement("div");
  tooltip.className = "chart-tooltip";
  tooltip.hidden = true;
  tooltip.setAttribute("role", "status");
  tooltip.setAttribute("aria-live", "off");
  wrapper.append(tooltip);
  let selected = Math.max(0, bins.findIndex(bin => candidate.original_energy >= bin.start && candidate.original_energy < bin.end));
  function show(index, clientX, clientY, keyboard = false) {
    selected = index;
    const bin = bins[index];
    const examples = bin.entries.slice(0, 3).map(item => `#${String(item.shuffle).padStart(3, "0")} ${energy(item.energy)}`).join(" · ");
    tooltip.setAttribute("aria-live", keyboard ? "polite" : "off");
    tooltip.innerHTML = `<strong>${bin.start} to ${bin.end} kcal/mol</strong><span>${bin.count} ${bin.count === 1 ? "control" : "controls"}</span>${examples ? `<small>${examples}${bin.count > 3 ? " …" : ""}</small>` : ""}`;
    tooltip.hidden = false;
    const bounds = wrapper.getBoundingClientRect();
    tooltip.style.left = `${Math.max(8, Math.min(bounds.width - 215, (clientX ?? bounds.left + bounds.width / 2) - bounds.left - 100))}px`;
    tooltip.style.top = `${Math.max(4, (clientY ?? bounds.top + 100) - bounds.top - 80)}px`;
    svg.querySelectorAll(".hist-bar").forEach((bar, position) => bar.classList.toggle("is-active", position === index));
  }
  function hide() {
    tooltip.hidden = true;
    svg.querySelectorAll(".hist-bar").forEach(bar => bar.classList.remove("is-active"));
  }
  svg.addEventListener("pointermove", event => {
    const bar = event.target.closest(".hist-bar");
    if (bar) show(Number(bar.dataset.bin), event.clientX, event.clientY);
    else hide();
  });
  svg.addEventListener("pointerleave", hide);
  svg.addEventListener("focus", () => show(selected, undefined, undefined, true));
  svg.addEventListener("blur", hide);
  svg.addEventListener("keydown", event => {
    if (event.key === "Escape") hide();
    if (["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) {
      event.preventDefault();
      if (event.key === "Home") selected = 0;
      else if (event.key === "End") selected = bins.length - 1;
      else selected = Math.max(0, Math.min(bins.length - 1, selected + (event.key === "ArrowRight" ? 1 : -1)));
      show(selected, undefined, undefined, true);
    }
  });
}

function renderCharts(data) {
  const edges = data.metadata.common_bin_edges;
  const axisMax = data.metadata.common_histogram_max_count;
  const featured = data.candidates.reduce((best, candidate) => candidate.tail_fraction < best.tail_fraction ? candidate : best);
  const hero = document.getElementById("hero-chart");
  hero.innerHTML = `<div class="chart-svg-wrap">${chartSVG(featured, edges, axisMax, "hero")}</div>`;
  attachChartTooltip(hero.querySelector(".chart-svg-wrap"), featured, edges);
  const grid = document.getElementById("distribution-grid");
  grid.innerHTML = data.candidates.map(candidate => `<article class="chart-card">
    <header class="chart-card-header"><div><p class="chart-kicker">${escapeHTML(geneLabel(candidate))} · ${candidate.length_nt} nt</p><h3><button class="candidate-link" data-select-candidate="${escapeHTML(candidate.id)}">${escapeHTML(candidate.id)}</button></h3></div><button class="chart-download" data-download-chart="${escapeHTML(candidate.id)}" aria-label="Download chart for ${escapeHTML(candidate.id)}">SVG <span aria-hidden="true">↗</span></button></header>
    <div class="chart-svg-wrap">${chartSVG(candidate, edges, axisMax, "panel")}</div>
    <footer class="chart-footer"><span><strong>${candidate.controls_le_original}/100</strong> equally or more favorable</span><span>Tail fraction <strong>${fraction(candidate.tail_fraction)}</strong></span></footer>
    ${candidate.related_group === "shared_626nt" ? '<p class="chart-meta">Related candidate · shared 626-nt block</p>' : '<p class="chart-meta">LINC00869 · smallest energy tail fraction</p>'}
  </article>`).join("") + `<article class="chart-guide"><p class="eyebrow">Reading the comparison</p><h3>Energy and control-relative evidence</h3><p>The four related candidates have the most negative original energies. LINC00869 has the smallest fraction of controls that perform as well.</p><div class="formula-note"><span>Exploratory energy tail fraction</span><strong>(equally or more favorable controls + 1) / 101</strong></div><p>All five fractions exceed 0.05. Every panel uses the same bins and axes.</p></article>`;
  grid.querySelectorAll(".chart-card").forEach((card, index) => attachChartTooltip(card.querySelector(".chart-svg-wrap"), data.candidates[index], edges));
  grid.addEventListener("click", event => {
    const download = event.target.closest("[data-download-chart]");
    if (download) {
      const candidate = data.candidates.find(item => item.id === download.dataset.downloadChart);
      downloadFile(`${candidate.id}-control-energies.svg`, chartSVG(candidate, edges, axisMax, "export"), "image/svg+xml");
    }
    const choose = event.target.closest("[data-select-candidate]");
    if (choose) {
      const select = document.getElementById("candidate-select");
      select.value = choose.dataset.selectCandidate;
      select.dispatchEvent(new Event("change"));
      document.getElementById("candidates").scrollIntoView({ behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "instant" : "smooth" });
      select.focus({ preventScroll: true });
    }
  });
  const downloadAll = document.getElementById("download-charts");
  downloadAll.addEventListener("click", () => {
    const panels = data.candidates.map((candidate, index) => {
      const x = 30 + (index % 2) * 550;
      const y = 100 + Math.floor(index / 2) * 335;
      const chart = chartSVG(candidate, edges, axisMax, `poster-${index}`).replace('<svg class="histogram"', `<svg x="${x}" y="${y + 28}" width="520" height="280" class="histogram"`);
      return `<text x="${x + 20}" y="${y + 17}" font-family="Arial, sans-serif" font-size="16" fill="${COLORS.ink}">${escapeHTML(candidate.id)} · ${escapeHTML(geneLabel(candidate))}</text>${chart}`;
    }).join("");
    const poster = `<svg xmlns="http://www.w3.org/2000/svg" width="1130" height="1145" viewBox="0 0 1130 1145"><rect width="1130" height="1145" fill="${COLORS.paper}"/><text x="50" y="42" font-family="Georgia, serif" font-size="28" fill="${COLORS.ink}">Human lncRNA–HPV16 interaction energies</text><text x="50" y="72" font-family="Arial, sans-serif" font-size="14" fill="${COLORS.muted}">Exploratory analysis · 100 actual shuffled controls per transcript · common bins and axes</text>${panels}<text x="50" y="1124" font-family="Arial, sans-serif" font-size="13" fill="${COLORS.muted}">IntaRNA 3.4.1 / ViennaRNA 2.7.0 · Coral line marks original energy · No experimental binding established</text></svg>`;
    downloadFile("BINF-five-control-distributions.svg", poster, "image/svg+xml");
  });
}

function renderTable(data) {
  const body = document.getElementById("comparison-body");
  let sortKey = "rank", descending = false;
  const sortButtons = [...document.querySelectorAll("[data-sort]")];
  function update() {
    const candidates = [...data.candidates].sort((a, b) => {
      const key = sortKey === "length" ? "length_nt" : sortKey;
      const comparison = typeof a[key] === "string" ? a[key].localeCompare(b[key]) : a[key] - b[key];
      return (descending ? -comparison : comparison) || a.rank - b.rank;
    });
    body.innerHTML = candidates.map(candidate => `<tr>
      <th scope="row"><a href="#candidates" data-table-candidate="${escapeHTML(candidate.id)}">${escapeHTML(candidate.id)}</a>${candidate.related_group === "shared_626nt" ? '<small>Shared 626-nt block</small>' : ''}</th>
      <td>${escapeHTML(geneLabel(candidate))}<small>${number.format(candidate.length_nt)} nt</small></td>
      <td class="numeric original-value">${energy(candidate.original_energy)}</td>
      <td class="numeric">${Number(candidate.shuffled_mean).toFixed(4)}</td>
      <td class="numeric">${Number(candidate.shuffled_median).toFixed(3)}</td>
      <td class="numeric">${candidate.controls_le_original}<span class="muted"> / 100</span></td>
      <td class="numeric tail-value">${fraction(candidate.tail_fraction)}</td>
    </tr>`).join("");
    sortButtons.forEach(button => {
      const active = button.dataset.sort === sortKey;
      button.closest("th").setAttribute("aria-sort", active ? (descending ? "descending" : "ascending") : "none");
      button.classList.toggle("is-sorted", active);
    });
  }
  sortButtons.forEach(button => button.addEventListener("click", () => {
    if (sortKey === button.dataset.sort) descending = !descending;
    else { sortKey = button.dataset.sort; descending = false; }
    update();
  }));
  body.addEventListener("click", event => {
    const link = event.target.closest("[data-table-candidate]");
    if (link) {
      const select = document.getElementById("candidate-select");
      select.value = link.dataset.tableCandidate;
      select.dispatchEvent(new Event("change"));
    }
  });
  update();
}

function coordinateTrack(start1, end1, length, label) {
  return `<div class="span-track" role="img" aria-label="${escapeHTML(label)} bases ${start1} through ${end1} of ${length}, one-based inclusive"><span class="span-fill" style="left:${100 * (start1 - 1) / length}%;width:${100 * (end1 - start1 + 1) / length}%"></span></div><div class="track-endpoints"><span>1</span><span>${number.format(length)} nt</span></div>`;
}

function renderCandidateDetails(data) {
  const select = document.getElementById("candidate-select");
  select.innerHTML = data.candidates.map(candidate => `<option value="${escapeHTML(candidate.id)}">${escapeHTML(candidate.id)} · ${escapeHTML(geneLabel(candidate))}</option>`).join("");
  function update() {
    const candidate = data.candidates.find(item => item.id === select.value);
    const [humanStart, humanEnd] = candidate.predicted_human;
    const [viralStart, viralEnd] = candidate.predicted_viral;
    const related = candidate.related_group === "shared_626nt";
    document.getElementById("candidate-detail").innerHTML = `<div class="candidate-intro"><p class="eyebrow">${escapeHTML(geneLabel(candidate))}</p><h3>${escapeHTML(candidate.id)}</h3><p>${number.format(candidate.length_nt)} nt · ${escapeHTML(candidate.gene_id)} · GENCODE v50</p></div>
      <div class="candidate-facts"><div class="fact"><span>Exact-match rows</span><strong>${candidate.step4.real_count}</strong></div><div class="fact"><span>Overlap clusters</span><strong>${candidate.overlap_clusters}</strong></div><div class="fact"><span>Nonoverlapping 10-mers</span><strong>${candidate.max_nonoverlapping_sites}</strong></div><div class="fact"><span>Step 4 BH-adjusted value</span><strong>${fraction(candidate.step4.bh_q_value)}</strong></div></div>
      <div class="detail-grid"><section class="coordinate-block"><h4>Original predicted regions</h4><p class="coordinate-convention">One-based, inclusive coordinates along the full RNAs.</p><div class="span-label"><span>Human RNA</span><strong>${humanStart}–${humanEnd}</strong></div>${coordinateTrack(humanStart, humanEnd, candidate.length_nt, "Human RNA predicted span")}<div class="span-label"><span>HPV16 RNA</span><strong>${viralStart}–${viralEnd}</strong></div>${coordinateTrack(viralStart, viralEnd, data.metadata.target.length_nt, "Viral RNA predicted span")}<p class="coordinate-note">Opening costs · human ${energy(candidate.ED1)}, viral ${energy(candidate.ED2)} kcal/mol</p></section>
      <section class="coordinate-block"><h4>Where the matches accumulate</h4><p><strong>${candidate.dominant_repeat.match_rows} of ${candidate.step4.real_count}</strong> exact-match rows fall within a <strong>${candidate.dominant_repeat.length_nt}-nt</strong> CCU-rich region.</p><p>Overlapping windows and repeated motifs contribute to this count. The ${candidate.max_nonoverlapping_sites} nonoverlapping 10-mers describe site geometry and can still share repeated sequence.</p><div class="relationship-note"><strong>${related ? "Part of the shared-sequence group" : "LINC00869 repeat context"}</strong><p>${related ? "This candidate shares an identical 626-nt block with three other shortlisted transcripts. All exact matches lie within that block, and the four original predictions use identical span sequences. Their control distributions are retained separately." : "Its full sequence differs from the other candidates. A 38-nt repeat region is also present within the four related transcripts."}</p></div></section></div>`;
  }
  select.addEventListener("change", update);
  update();
}

function verifyDataShape(data) {
  if (data.candidates.length !== 5 || new Set(data.candidates.map(candidate => candidate.id)).size !== 5) throw new Error("Expected five distinct transcript records.");
  data.candidates.forEach(candidate => {
    if (candidate.energies.length !== 100) throw new Error("Each transcript needs 100 actual control energies.");
    histogram(candidate.energies, data.metadata.common_bin_edges);
    const count = candidate.energies.filter(value => value <= candidate.original_energy).length;
    if (count !== candidate.controls_le_original || Math.abs((count + 1) / 101 - candidate.tail_fraction) > 1e-10) throw new Error("The prepared energy statistics are inconsistent.");
  });
}

async function initialize() {
  try {
    const response = await fetch("data/website-data.json");
    if (!response.ok) throw new Error(`Data request failed (${response.status}).`);
    const data = await response.json();
    verifyDataShape(data);
    renderCharts(data);
    renderTable(data);
    renderCandidateDetails(data);
    document.getElementById("provenance-details").textContent = JSON.stringify({ software: data.metadata.software, sources: data.metadata.sources, definitions: data.definitions }, null, 2);
    const status = document.getElementById("site-status");
    if (status) { status.textContent = "Verified data loaded. Five candidates, 500 control energies."; status.classList.add("visually-hidden"); }
    document.body.dataset.ready = "true";
  } catch (error) {
    console.error(error);
    const status = document.getElementById("site-status") || document.createElement("p");
    status.id = "site-status";
    status.className = "data-error";
    status.setAttribute("role", "alert");
    status.textContent = `The verified results could not be loaded. Serve this folder with the local server described in README.md. ${error.message}`;
    if (!status.isConnected) document.querySelector("main").prepend(status);
    document.body.dataset.ready = "error";
  }
}

if (typeof module !== "undefined") module.exports = { histogram, chartSVG, verifyDataShape };
if (typeof document !== "undefined") initialize();
