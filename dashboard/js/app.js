// Shell: opens the release, routes on the hash, runs the top bar and the search box. Pages live in their own modules.
import { ensure, query, lit, openRelease, releaseList, openEntry, NoRelease } from "./db.js";
import { setOverview, state, esc, NETWORKS, networkHref } from "./format.js";
import * as landing from "./landing.js";
import * as overview from "./overview.js";
import * as term from "./term.js";
import * as gene from "./gene.js";
import * as phenotype from "./phenotype.js";
import * as family from "./family.js";
import * as species from "./species.js";
import * as browse from "./browse.js";
import * as search from "./search.js";
import * as candidates from "./candidates.js";
import * as downloads from "./downloads.js";

const ROUTES = [
  [/^#\/?$/, landing.render],
  [/^#\/release$/, overview.render],
  [/^#\/(?:disease|human-phenotype)\/(.+)$/, term.render],
  [/^#\/gene\/(.+)$/, gene.render],
  [/^#\/phenotype\/(.+)$/, phenotype.render],
  [/^#\/family\/(.+)$/, family.render],
  [/^#\/species\/(.+)$/, species.render],
  [/^#\/browse\/([\w-]+)$/, browse.render],
  [/^#\/candidates$/, candidates.render],
  [/^#\/downloads$/, downloads.render],
  [/^#\/search\/(.*)$/, search.render],
];

// Addresses that belong to one network: a disease or human phenotype page, and the list of them
const NETWORK_PAGE = /^#\/(?:browse\/)?(diseases?|human-phenotypes?)(?:\/|$)/;

const view = document.getElementById("view");

function reportError(err) {
  const box = document.getElementById("errors");
  box.hidden = false;
  box.textContent += (err && err.stack ? err.stack : String(err)) + "\n";
  console.error(err);
}
window.addEventListener("error", e => reportError(e.error || e.message));
window.addEventListener("unhandledrejection", e => reportError(e.reason));

async function navigate() {
  const hash = location.hash || "#/";

  // a link into the other network's release opens that release
  const m = hash.match(NETWORK_PAGE);
  const wanted = m && (m[1].startsWith("disease") ? "disease" : "phenotype");
  if (wanted && wanted !== state.net.key) { location.replace(networkHref(wanted, hash)); return; }
  markNav(hash);

  for (const [pattern, render] of ROUTES) {
    const m = hash.match(pattern);
    if (!m) continue;
    const params = m.slice(1).map(decodeURIComponent);
    view.innerHTML = `<div class="loading">Loading…</div>`;
    window.scrollTo(0, 0);
    try {
      await render(view, ...params);
    } catch (err) {
      view.innerHTML = `<div class="card"><h3>Something went wrong</h3><p class="muted">${esc(err.message || err)}</p></div>`;
      reportError(err);
    }
    return;
  }
  view.innerHTML = `<div class="card"><h3>Not found</h3><p class="muted">${esc(hash)}</p></div>`;
}

// ---- top bar: a row per network with its pages, in columns both rows share ------------------
// The open network's links stay in this page. The other's reload into that network's release of the same KG
// version, and are dimmed when there is none
const PAGES = net => [["Release", "#/release"], ["Candidates", "#/candidates"], [net.Terms, `#/browse/${net.browse}`], ["Genes", "#/browse/genes"],
                      ["Phenotypes", "#/browse/phenotypes"], ["Families", "#/browse/families"], ["Species", "#/browse/species"]];

function networkRows() {
  const version = openEntry().graph_version;
  document.getElementById("nets").innerHTML = Object.values(NETWORKS).map(net => {
    const exists = releaseList().some(r => r.graph_version === version && r.prediction_network === net.key);
    const href = hash => !exists ? "" : `href="${net.key === state.net.key ? hash : networkHref(net.key, hash)}"`;
    const title = exists ? `${net.humans} against the phenotypes of each model species` : `No ${net.key} release for Monarch KG ${version}`;
    return `<div class="net-row${exists ? "" : " off"}" data-net="${net.key}"><span class="net-label" title="${esc(title)}">${esc(net.title)}</span>` +
      PAGES(net).map(([label, hash]) => `<a ${href(hash)} data-hash="${hash}">${esc(label)}</a>`).join("") + `</div>`;
  }).join("");
}

// The open page's network row, and its section there (a gene's page is under Genes). The landing page and
// Downloads cover both networks
const SECTION = { gene: "genes", phenotype: "phenotypes", family: "families", species: "species", disease: "diseases", "human-phenotype": "human-phenotypes" };

function markNav(hash) {
  const parts = hash.split("/");
  const page = SECTION[parts[1]] ? `#/browse/${SECTION[parts[1]]}` : parts.slice(0, 3).join("/");
  const shared = /^#\/?(downloads)?$/.test(hash);
  document.querySelectorAll(".net-row").forEach(row => {
    const here = !shared && row.dataset.net === state.net.key;
    row.classList.toggle("here", here);
    row.querySelectorAll("a").forEach(a => a.classList.toggle("active", here && a.dataset.hash === page));
  });
  document.getElementById("nav-downloads").classList.toggle("active", hash === "#/downloads");
}

// ---- search box: typeahead over the search table ---------------------------------------------
function setupSearch() {
  const form = document.getElementById("search-form");
  const input = document.getElementById("search-input");
  const box = document.getElementById("search-suggest");
  let timer = null, selected = -1, items = [];

  const hide = () => { box.hidden = true; selected = -1; };
  const show = () => {
    box.innerHTML = items.map((r, i) => `
      <a href="${search.href(r)}" class="${i === selected ? "sel" : ""}">
        <span class="kind">${r.kind === "term" ? state.net.term : r.kind}</span><span>${esc(r.label)} <span class="id muted">${esc(r.id)}</span></span>
        <span class="detail">${esc(r.detail || "")}</span></a>`).join("");
    box.hidden = items.length === 0;
  };

  input.addEventListener("input", () => {
    clearTimeout(timer);
    const q = input.value.trim();
    if (q.length < 2) { hide(); return; }
    timer = setTimeout(async () => {
      await ensure("search");
      const like = lit("%" + q + "%"), prefix = lit(q + "%");
      items = await query(`
        SELECT kind, id, label, detail FROM search
        WHERE label ILIKE ${like} OR id ILIKE ${prefix}
        ORDER BY (id ILIKE ${prefix}) DESC, (label ILIKE ${prefix}) DESC, weight DESC, label
        LIMIT 10`);
      selected = -1;
      show();
    }, 150);
  });
  input.addEventListener("keydown", e => {
    if (box.hidden) return;
    if (e.key === "ArrowDown") { selected = Math.min(selected + 1, items.length - 1); show(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { selected = Math.max(selected - 1, 0); show(); e.preventDefault(); }
    else if (e.key === "Escape") hide();
  });
  form.addEventListener("submit", e => {
    e.preventDefault();
    if (selected >= 0 && items[selected]) location.hash = search.href(items[selected]);
    else if (input.value.trim()) location.hash = `#/search/${encodeURIComponent(input.value.trim())}`;
    hide();
    input.blur();
  });
  document.addEventListener("click", e => { if (!form.contains(e.target)) hide(); });
}

// Nothing to open: say so, and list the releases there are
function noRelease(message) {
  const list = releaseList();
  const href = r => { const p = new URLSearchParams(location.search); p.set("release", r.graph_version); p.set("network", r.prediction_network); return `${location.pathname}?${p}#/`; };
  view.innerHTML = `<div class="card"><h3>${esc(message)}</h3>
    ${list.length ? `<p class="muted">The releases there are:</p><div class="name-grid wide">${list.map(r =>
        `<a href="${href(r)}">Monarch KG ${esc(r.graph_version)} · ${esc(NETWORKS[r.prediction_network]?.title.toLowerCase() || r.prediction_network)}</a>`).join("")}</div>`
      : `<p class="muted">Releases show up here as soon as one is published. How they are made: <a href="https://github.com/monarch-initiative/monarch-phenologs">monarch-phenologs</a>.</p>`}
  </div>`;
}

async function boot() {
  let ov;
  try {
    ov = await openRelease();
  } catch (err) {
    if (err instanceof NoRelease) { noRelease(err.message); return; }
    throw err;
  }
  setOverview(ov);
  const r = ov.release, net = state.net;
  document.getElementById("brand-release").textContent = `kg ${r.graph_version}`;
  document.getElementById("nav-downloads").hidden = false;
  networkRows();
  document.getElementById("search-input").placeholder = `Search a ${net.term}, gene, phenotype, protein family or species…`;
  document.getElementById("footer-release").textContent =
    `Monarch KG ${r.graph_version}${r.phenio_version ? " · phenio " + r.phenio_version : ""} · ${net.title.toLowerCase()} · FDR ${r.fdr} · ` +
    `${r.rank_metric === "hg" ? "hypergeometric" : "naive Bayes"} ranking over each ${net.term}'s ${Number(r.kneighbs).toLocaleString("en-US")} nearest phenologs · built ${r.built}`;
  setupSearch();
  window.addEventListener("hashchange", navigate);
  await navigate();
  ensure("search");           // warm the index behind the search box
}

boot().catch(reportError);
