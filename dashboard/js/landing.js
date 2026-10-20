// The landing page: one Monarch KG version's two releases side by side (disease phenologs and phenotype phenologs),
// their model species in both, and every release there is
import { releaseList, openEntry, overviewOf } from "./db.js";
import { fmt, esc, NETWORKS, networkHref, route } from "./format.js";

const KEYS = ["disease", "phenotype"];

export async function render(el) {
  const version = openEntry().graph_version;
  const entries = KEYS.map(k => releaseList().find(r => r.graph_version === version && r.prediction_network === k));
  const ovs = await Promise.all(entries.map(e => e ? overviewOf(e) : null));
  const first = ovs.find(Boolean);
  const r = first.release;
  const model = first.species.filter(s => !s.is_base);

  el.innerHTML = `
    <section class="card">
      <div class="kicker">Monarch KG release</div>
      <h1 class="landing-title">Phenologs from the Monarch KG <span class="mono" style="color:var(--accent)">${esc(version)}</span></h1>
      <p class="lede">Orthologous phenotypes between human and ${model.length} model species, found two ways: from human diseases
        (Mondo) and from human phenotypes (HPO), each compared with every model species' phenotypes through the PANTHER protein
        families of their genes. Each disease or human phenotype then ranks the human genes its phenologs point to as candidates.</p>
      <div class="meta"><span>graph <b>${esc(version)}</b></span>${r.phenio_version ? `<span>phenio <b>${esc(r.phenio_version)}</b></span>` : ""}
        ${r.graph_url ? `<span><a href="${esc(r.graph_url)}" target="_blank" rel="noopener">the graph ↗</a></span>` : ""}</div>
    </section>

    <section class="section grid cols-2">${KEYS.map((k, i) => networkCard(k, ovs[i])).join("")}</section>

    <section class="section card">
      <div class="card-head"><h2>Model species</h2><span class="sub">phenologs with each, in each network</span></div>
      <div id="species"></div>
    </section>

    <section class="section card">
      <div class="card-head"><h2>Releases</h2><span class="sub">one per Monarch KG version and network</span></div>
      <div id="releases"></div>
    </section>`;

  speciesTable(el.querySelector("#species"), model, ovs);
  releasesTable(el.querySelector("#releases"), version);
}

// A network's headline numbers, with the way into its release
function networkCard(key, ov) {
  const net = NETWORKS[key];
  if (!ov) {
    return `<div class="card net-card"><div class="card-head"><h2>${esc(net.title)}</h2></div>
      <div class="empty">No ${esc(net.key)} release for this KG version</div></div>`;
  }
  const t = ov.tiles, r = ov.release;
  return `
    <div class="card net-card">
      <div class="card-head"><h2>${esc(net.title)}</h2><span class="sub">${esc(net.humans.toLowerCase())} × model species phenotypes</span></div>
      <div class="tiles net-tiles">
        ${tile(t.terms, `${net.terms} with phenologs`, "accent")}
        ${tile(t.phenologs, "significant phenologs", "accent")}
        ${tile(t.model_phenotypes, `model phenotypes with ${net.terms}`)}
        ${tile(t.scored_pairs, `ranked ${net.pair} pairs`)}
        ${tile(t.known_recovered, `known ${net.pair} pairs recovered (${fmt.pct(t.known_recovered / t.known_in_reach)} within reach)`)}
        ${tile(t.candidate_genes, `genes with one or more candidate ${net.term} associations`)}
      </div>
      <div class="net-foot">
        <span class="faint">fdr ${r.fdr} · k ${fmt.int(r.kneighbs)} · ${esc(r.rank_metric)}${r.top_candidates ? ` · best ${fmt.int(r.top_candidates)} candidates kept per ${esc(net.term)}` : ""}
          · built ${esc(r.built)}${r.code_version ? ` · code ${esc(r.code_version)}` : ""}</span>
        <span class="go"><a href="${networkHref(key, "#/downloads")}">Downloads</a><a class="btn" href="${networkHref(key, "#/release")}">Open ${esc(net.title.toLowerCase())} →</a></span>
      </div>
    </div>`;
}

function tile(n, label, cls = "") {
  return `<div class="tile ${cls}"><div class="n">${fmt.int(n)}</div><div class="l">${esc(label)}</div></div>`;
}

// One row per model species: its phenologs and terms in each network, a hairline bar under each count (scaled within its column)
function speciesTable(el, model, ovs) {
  const cols = KEYS.flatMap((k, i) => ovs[i] ? [[i, "n_phenologs", `${NETWORKS[k].title}`], [i, "n_terms", NETWORKS[k].Terms]] : []);
  const value = (i, key, sp) => ovs[i].species.find(s => s.species === sp)?.[key] ?? 0;
  const max = cols.map(([i, key]) => Math.max(1, ...model.map(s => value(i, key, s.species))));
  el.innerHTML = `<div class="tbl-wrap"><table class="tbl inputs wrap-head"><thead><tr><th>Species</th>
      ${cols.map(([, , label]) => `<th class="num">${esc(label)}</th>`).join("")}</tr></thead><tbody>
    ${model.map(s => `<tr><td><a href="${route.species(s.species)}"><span class="dot" style="background:${s.color}"></span> <i>${esc(s.binomial)}</i></a>
        <div class="faint">${esc(s.common)}</div></td>
      ${cols.map(([i, key], c) => { const v = value(i, key, s.species);
        return `<td class="num">${fmt.int(v)}<span class="hair" style="width:${(100 * v / max[c]).toFixed(1)}%; background:${s.color}"></span></td>`; }).join("")}
    </tr>`).join("")}</tbody></table></div>`;
}

// Every KG version with a release, newest first, and its networks
function releasesTable(el, current) {
  const byVersion = d3.group(releaseList(), r => r.graph_version);
  const href = (version, network, hash) => {
    const p = new URLSearchParams(location.search);
    p.set("release", version);
    if (network) p.set("network", network);
    return `${location.pathname}?${p}${hash}`;
  };
  el.innerHTML = `<div class="tbl-wrap"><table class="tbl"><thead><tr><th>Monarch KG</th><th>Phenio</th><th>Releases</th><th>Built</th><th>Code</th></tr></thead><tbody>
    ${[...byVersion.entries()].map(([version, rs]) => `<tr class="${version === current ? "current-row" : ""}">
      <td class="mono"><a href="${href(version, null, "#/")}">${esc(version)}</a></td>
      <td class="mono">${esc(rs[0].phenio_version || "")}</td>
      <td>${KEYS.filter(k => rs.some(r => r.prediction_network === k)).map(k => `<a href="${href(version, k, "#/release")}">${esc(NETWORKS[k].title.toLowerCase())}</a>`).join(" · ")}</td>
      <td>${esc([...new Set(rs.map(r => r.built))].join(", "))}</td>
      <td class="mono">${esc([...new Set(rs.map(r => r.code_version || ""))].join(", "))}</td></tr>`).join("")}
    </tbody></table></div>`;
}
