// The open release's files: what each holds, its rows and size, and how to read a table straight from where it is
import { getJSON, fileURL, ROOT, BUCKET } from "./db.js";
import { state, fmt, esc } from "./format.js";

export async function render(el) {
  const net = state.net, r = state.overview.release;
  const files = await getJSON(fileURL("manifest.json"));
  const bytes = files.reduce((n, f) => n + f.bytes, 0);
  const example = fileURL("candidates.parquet");

  el.innerHTML = `
    <div class="page-head"><div class="titles"><div class="kicker">Downloads</div>
      <h1>${esc(net.title)} · Monarch KG <span class="mono">${esc(r.graph_version)}</span></h1>
      <div class="muted">Every table the site shows is a parquet file of this release. <a href="${fileURL("README.md")}">README.md</a>
        describes each table and column; <a href="${fileURL("manifest.json")}">manifest.json</a> lists the files with their sha256.</div></div></div>
    <div class="card"><div class="card-head"><h2>Files</h2><span class="sub">${fmt.int(files.length)} files · ${fmt.dec(bytes / 1e6)} MB</span></div>
      <div class="tbl-wrap"><table class="tbl"><thead><tr><th>File</th><th>What it holds</th><th class="num">Rows</th><th class="num">Size</th></tr></thead><tbody>
      ${files.map(f => `<tr><td class="mono nowrap"><a href="${fileURL(f.path)}" download>${esc(f.path)}</a></td><td>${esc(f.description)}</td>
        <td class="num">${f.rows == null ? "" : fmt.int(f.rows)}</td><td class="num nowrap">${size(f.bytes)}</td></tr>`).join("")}
      </tbody></table></div></div>
    <div class="card section"><div class="card-head"><h2>Reading a table where it is</h2><span class="sub">DuckDB, Python and R, straight from its address</span></div>
      <pre class="code">duckdb -c "SELECT * FROM '${esc(example)}' WHERE NOT known ORDER BY global_rank LIMIT 10"
python -c "import pandas as pd; print(pd.read_parquet('${esc(example)}').head())"
Rscript -e 'download.file("${esc(example)}", f &lt;- tempfile(), mode = "wb"); print(head(arrow::read_parquet(f)))'</pre>
      ${ROOT === BUCKET ? `<p class="muted" style="margin-bottom:0">For scripts that should always read the newest ${esc(net.term)} release:
        <span class="mono">${esc(BUCKET + "latest/" + net.key + "/")}</span></p>` : ""}
    </div>`;
}

function size(bytes) {
  return bytes < 1e6 ? `${fmt.int(Math.max(1, Math.round(bytes / 1e3)))} kB` : `${fmt.dec(bytes / 1e6)} MB`;
}
