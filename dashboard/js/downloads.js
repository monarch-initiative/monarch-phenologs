// The downloads of a Monarch KG version: each network's phenologs and gene candidates, the networks' READMEs, and a
// command per file that downloads it as .tsv
import { getJSON, fileURL, fileOf, releaseList, openEntry, ROOT, BUCKET } from "./db.js";
import { fmt, esc, NETWORKS } from "./format.js";

export async function render(el) {
  const version = openEntry().graph_version;
  const entries = Object.keys(NETWORKS).map(k => releaseList().find(r => r.graph_version === version && r.prediction_network === k)).filter(Boolean);
  const manifests = await Promise.all(entries.map(e => getJSON(fileURL("manifest.json", e))));

  // A file of a network's release, with what its manifest says of it
  const file = (i, path, name = path) => {
    const f = manifests[i].find(m => m.path === path) || {};
    return { name, url: fileURL(path, entries[i]), description: f.description, rows: f.rows, bytes: f.bytes };
  };
  const tables = ["phenologs", "candidates"].flatMap(t => entries.map((e, i) => file(i, fileOf(t, e))));
  const readmes = entries.map((e, i) => ({ ...file(i, "README.md", `${e.prediction_network}/README.md`),
                                           description: `What every table and column of the ${NETWORKS[e.prediction_network].title.toLowerCase()} holds` }));
  const files = [...tables, ...readmes];
  const tsv = f => `duckdb -c "COPY (SELECT * FROM '${f.url}') TO '${f.name.replace(/\.parquet$/, ".tsv")}' (DELIMITER '\\t')"`;

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Downloads</div>
      <h1>Phenologs · Monarch KG <span class="mono">${esc(version)}</span></h1>
      <div class="muted">The phenologs and gene candidates of each network, as parquet. Each network's README describes every table and column.</div></div></div>
    <div class="card"><div class="card-head"><h2>Files</h2><span class="sub">${fmt.int(files.length)} files · ${fmt.dec(d3.sum(files, f => f.bytes) / 1e6)} MB</span></div>
      <div class="tbl-wrap"><table class="tbl"><thead><tr><th class="fixed">File</th><th class="fixed">What it holds</th><th class="num fixed">Rows</th><th class="num fixed">Size</th></tr></thead><tbody>
      ${files.map(f => `<tr><td class="mono nowrap"><a href="${esc(f.url)}">${esc(f.name)}</a></td><td>${esc(f.description)}</td>
        <td class="num">${f.rows == null ? "" : fmt.int(f.rows)}</td><td class="num nowrap">${f.bytes == null ? "" : size(f.bytes)}</td></tr>`).join("")}
      </tbody></table></div></div>
    <div class="card section"><div class="card-head"><h2>Download and convert to .tsv</h2>
        <span class="sub">one command per file, with <a href="https://duckdb.org/docs/installation/" target="_blank" rel="noopener">DuckDB</a></span></div>
      <div class="cmds">${tables.map(f => `<div class="cmd"><pre class="code">${esc(tsv(f))}</pre>${navigator.clipboard ? `<button class="copy">Copy</button>` : ""}</div>`).join("")}</div>
      ${ROOT === BUCKET ? `<p class="muted" style="margin:10px 0 0">For scripts that should always get the newest release, put
        <span class="mono">latest/</span> in place of <span class="mono">${esc(version)}/</span>.</p>` : ""}
    </div>`;

  el.querySelectorAll("button.copy").forEach(b => b.addEventListener("click", () =>
    navigator.clipboard.writeText(b.previousElementSibling.textContent).then(() => {
      b.textContent = "Copied";
      setTimeout(() => { b.textContent = "Copy"; }, 1500);
    }, () => {})));
}

function size(bytes) {
  return bytes < 1e6 ? `${fmt.int(Math.max(1, Math.round(bytes / 1e3)))} kB` : `${fmt.dec(bytes / 1e6)} MB`;
}
