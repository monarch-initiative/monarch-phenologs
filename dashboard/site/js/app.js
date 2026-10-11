// Shell: loads overview.json, routes on the hash, runs the search box. Views live in ./views.
import { ensure, query, lit } from "./db.js";
import { setOverview, state, esc } from "./format.js";
import * as overview from "./views/overview.js";
import * as disease from "./views/disease.js";
import * as gene from "./views/gene.js";
import * as phenotype from "./views/phenotype.js";
import * as family from "./views/family.js";
import * as species from "./views/species.js";
import * as browse from "./views/browse.js";
import * as search from "./views/search.js";

const ROUTES = [
  [/^#\/?$/, overview.render, ""],
  [/^#\/disease\/(.+)$/, disease.render],
  [/^#\/gene\/(.+)$/, gene.render],
  [/^#\/phenotype\/(.+)$/, phenotype.render],
  [/^#\/family\/(.+)$/, family.render],
  [/^#\/species\/(.+)$/, species.render],
  [/^#\/browse\/(\w+)$/, browse.render],
  [/^#\/search\/(.*)$/, search.render],
];

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

// ---- search box: typeahead over the search table ---------------------------------------------
function setupSearch() {
  const form = document.getElementById("search-form");
  const input = document.getElementById("search-input");
  const box = document.getElementById("search-suggest");
  let timer = null, selected = -1, items = [];

  const hide = () => { box.hidden = true; selected = -1; };
  const show = () => {
    box.innerHTML = items.map((r, i) => `
      <a href="#/${r.kind}/${encodeURIComponent(r.id)}" class="${i === selected ? "sel" : ""}">
        <span class="kind">${r.kind}</span><span>${esc(r.label)} <span class="id muted">${esc(r.id)}</span></span>
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
    if (selected >= 0 && items[selected]) location.hash = `#/${items[selected].kind}/${encodeURIComponent(items[selected].id)}`;
    else if (input.value.trim()) location.hash = `#/search/${encodeURIComponent(input.value.trim())}`;
    hide();
    input.blur();
  });
  document.addEventListener("click", e => { if (!form.contains(e.target)) hide(); });
}

async function boot() {
  const ov = await (await fetch("data/overview.json")).json();
  setOverview(ov);
  const r = ov.release;
  document.getElementById("brand-release").textContent = `kg ${r.graph_version}`;
  document.getElementById("footer-release").textContent =
    `Monarch KG ${r.graph_version}${r.phenio_version ? " · phenio " + r.phenio_version : ""} · FDR ${r.fdr} · k=${r.kneighbs} ${r.rank_metric} · built ${r.built}`;
  setupSearch();
  window.addEventListener("hashchange", navigate);
  await navigate();
  ensure("search");           // warm the index behind the search box
}

boot().catch(reportError);
