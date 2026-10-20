// Shell: opens the release, routes on the hash, runs the network switch and the search box. Pages live in their own modules.
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
  networkSwitch(hash);

  for (const [pattern, render] of ROUTES) {
    const m = hash.match(pattern);
    if (!m) continue;
    const params = m.slice(1).map(decodeURIComponent);
    document.querySelectorAll(".nav a").forEach(a => a.classList.toggle("active", a.getAttribute("href") === hash.split("/").slice(0, 3).join("/")));
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

// The disease / phenotype switch: this page in the other network's release of the same KG version. A disease or
// human phenotype page has no counterpart there, so it leads to that release's page
function networkSwitch(hash) {
  const version = openEntry().graph_version;
  document.querySelectorAll("#net-switch a").forEach(a => {
    const net = NETWORKS[a.dataset.net];
    const exists = releaseList().some(r => r.graph_version === version && r.prediction_network === net.key);
    const here = net.key === state.net.key;
    a.classList.toggle("active", here);
    a.classList.toggle("off", !exists);
    a.title = exists ? `${net.title}: ${net.humans.toLowerCase()} against the phenotypes of each model species` : `No ${net.key} release for Monarch KG ${version}`;
    const target = /^#\/(disease|human-phenotype)\//.test(hash) ? "#/release" : hash.match(NETWORK_PAGE) ? `#/browse/${net.browse}` : hash;
    if (here || !exists) a.removeAttribute("href"); else a.href = networkHref(net.key, target);
  });
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
  const terms = document.getElementById("nav-terms");
  terms.href = `#/browse/${net.browse}`;
  terms.textContent = net.Terms;
  document.getElementById("search-input").placeholder = `Search a ${net.term}, gene, phenotype, protein family or species…`;
  document.getElementById("net-switch").hidden = false;
  document.getElementById("footer-release").textContent =
    `Monarch KG ${r.graph_version}${r.phenio_version ? " · phenio " + r.phenio_version : ""} · ${net.title.toLowerCase()} · FDR ${r.fdr} · ` +
    `${r.rank_metric === "hg" ? "hypergeometric" : "naive Bayes"} ranking over each ${net.term}'s ${Number(r.kneighbs).toLocaleString("en-US")} nearest phenologs · built ${r.built}`;
  setupSearch();
  window.addEventListener("hashchange", navigate);
  await navigate();
  ensure("search");           // warm the index behind the search box
}

boot().catch(reportError);
