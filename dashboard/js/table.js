// Sortable, pageable, filterable tables.
//   table(el, rows, columns, options)            rows already in memory
//   pagedTable(el, columns, fetchPage, options)  rows too many to hold: every page, sort and filter is a query.
//                                                fetchPage({ sortKey, sortDir, filter, page, pageSize }) -> { rows, total }
//   columns: [{ key, label, num, format(value, row) -> html, sortValue(row), sortDir: first click, width, help: header tooltip, sortable }]
//   options: { pageSize, sortKey, sortDir, search: [keys to match], placeholder, tools: html, empty, wrapHead, onDraw(el) }
import { esc, fmt } from "./format.js";

export function table(el, rows, columns, options = {}) {
  const st = tableState(options);

  function sorted() {
    let out = rows;
    if (st.filter && options.search) {
      const q = st.filter.toLowerCase();
      out = out.filter(r => options.search.some(k => String(r[k] ?? "").toLowerCase().includes(q)));
    }
    if (st.sortKey) {
      const col = columns.find(c => c.key === st.sortKey) || {};
      const val = col.sortValue || (r => r[st.sortKey]);
      const dir = st.sortDir === "asc" ? 1 : -1;
      out = [...out].sort((a, b) => {
        const x = val(a), y = val(b);
        if (x == null && y == null) return 0;
        if (x == null) return 1;
        if (y == null) return -1;
        return (x < y ? -1 : x > y ? 1 : 0) * dir;
      });
    }
    return out;
  }

  function draw() {
    const all = sorted();
    st.page = Math.min(st.page, Math.max(0, Math.ceil(all.length / st.pageSize) - 1));
    render(el, st, columns, all.slice(st.page * st.pageSize, (st.page + 1) * st.pageSize), all.length, options, draw);
  }

  draw();
  return { redraw: draw, setRows(next) { rows = next; st.page = 0; draw(); } };
}

export function pagedTable(el, columns, fetchPage, options = {}) {
  const st = tableState(options);
  let latest = 0;

  async function draw() {
    const mine = ++latest;
    const { rows, total } = await fetchPage({ ...st });
    if (mine === latest) render(el, st, columns, rows, total, { debounce: 250, ...options }, draw);    // a newer request wins
  }

  draw();
  return { reload() { st.page = 0; return draw(); } };
}

function tableState(options) {
  return { page: 0, pageSize: options.pageSize || 25, sortKey: options.sortKey || null, sortDir: options.sortDir || "desc", filter: "", timer: null };
}

// One page of rows, with the header, pager and filter box wired back to redraw()
function render(el, st, columns, slice, total, options, redraw) {
  const pages = Math.max(1, Math.ceil(total / st.pageSize));
  const typing = el.contains(document.activeElement) && document.activeElement.matches(".tbl-tools input[type=search]");

  const head = columns.map(c => `
    <th class="${c.num ? "num" : ""} ${st.sortKey === c.key ? "sorted" : ""} ${c.sortable === false ? "fixed" : ""}" ${c.sortable === false ? "" : `data-key="${c.key}"`}
        ${c.width ? `style="width:${c.width}"` : ""} ${c.help ? `title="${esc(c.help)}"` : ""}>
      ${esc(c.label)}${st.sortKey === c.key ? (st.sortDir === "asc" ? " ▲" : " ▼") : ""}</th>`).join("");
  const body = slice.map(r => "<tr>" + columns.map(c => {
    const v = r[c.key];
    const html = c.format ? c.format(v, r) : c.num ? fmt.int(v) : esc(v ?? "");
    return `<td class="${c.num ? "num" : ""} ${c.ellipsis ? "ellipsis" : ""}" ${c.ellipsis ? `title="${esc(v ?? "")}"` : ""}>${html}</td>`;
  }).join("") + "</tr>").join("");

  el.innerHTML = `
    <div class="tbl-tools">
      ${options.search ? `<input type="search" placeholder="${esc(options.placeholder || "Filter…")}" value="${esc(st.filter)}">` : ""}
      ${options.tools || ""}
      <span class="muted" style="margin-left:auto">${total.toLocaleString("en-US")} rows</span>
    </div>
    <div class="tbl-wrap"><table class="tbl ${options.wrapHead ? "wrap-head" : ""}"><thead><tr>${head}</tr></thead>
    <tbody>${body || `<tr><td colspan="${columns.length}" class="empty">${esc(options.empty || "Nothing to show")}</td></tr>`}</tbody></table></div>
    ${pages > 1 ? `<div class="tbl-foot">
      <button data-page="prev" ${st.page === 0 ? "disabled" : ""}>‹ Prev</button>
      <span>Page ${(st.page + 1).toLocaleString("en-US")} of ${pages.toLocaleString("en-US")}</span>
      <button data-page="next" ${st.page >= pages - 1 ? "disabled" : ""}>Next ›</button></div>` : ""}`;

  el.querySelectorAll("th[data-key]").forEach(th => th.addEventListener("click", () => {
    const key = th.dataset.key;
    if (st.sortKey === key) st.sortDir = st.sortDir === "asc" ? "desc" : "asc";
    else { const col = columns.find(c => c.key === key) || {}; st.sortKey = key; st.sortDir = col.sortDir || (col.num ? "desc" : "asc"); }
    st.page = 0;
    redraw();
  }));
  el.querySelectorAll("button[data-page]").forEach(b => b.addEventListener("click", () => {
    st.page += b.dataset.page === "next" ? 1 : -1;
    redraw();
  }));
  const search = el.querySelector(".tbl-tools input[type=search]");
  if (search) {
    search.addEventListener("input", () => {
      st.filter = search.value;
      st.page = 0;
      clearTimeout(st.timer);
      if (options.debounce) st.timer = setTimeout(redraw, options.debounce); else redraw();
    });
    if (typing) { search.focus(); search.setSelectionRange(search.value.length, search.value.length); }
  }
  if (options.onDraw) options.onDraw(el);
}
