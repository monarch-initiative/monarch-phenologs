// A sortable, pageable, filterable table rendered from rows already in memory.
//   table(el, rows, columns, options)
//   columns: [{ key, label, num, format(value, row) -> html, sortValue(row), width, help: header tooltip }]
//   options: { pageSize, sortKey, sortDir, search: [keys to match], placeholder, tools: html, empty }
import { esc } from "../format.js";

export function table(el, rows, columns, options = {}) {
  const st = {
    page: 0,
    pageSize: options.pageSize || 25,
    sortKey: options.sortKey || null,
    sortDir: options.sortDir || "desc",
    filter: "",
  };

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
    const pages = Math.max(1, Math.ceil(all.length / st.pageSize));
    st.page = Math.min(st.page, pages - 1);
    const slice = all.slice(st.page * st.pageSize, (st.page + 1) * st.pageSize);

    const head = columns.map(c => `
      <th class="${c.num ? "num" : ""} ${st.sortKey === c.key ? "sorted" : ""}" data-key="${c.key}" ${c.width ? `style="width:${c.width}"` : ""} ${c.help ? `title="${esc(c.help)}"` : ""}>
        ${esc(c.label)}${st.sortKey === c.key ? (st.sortDir === "asc" ? " ▲" : " ▼") : ""}</th>`).join("");
    const body = slice.map(r => "<tr>" + columns.map(c => {
      const v = r[c.key];
      const html = c.format ? c.format(v, r) : esc(v ?? "");
      return `<td class="${c.num ? "num" : ""} ${c.ellipsis ? "ellipsis" : ""}" ${c.ellipsis ? `title="${esc(v ?? "")}"` : ""}>${html}</td>`;
    }).join("") + "</tr>").join("");

    el.innerHTML = `
      <div class="tbl-tools">
        ${options.search ? `<input type="search" placeholder="${esc(options.placeholder || "Filter…")}" value="${esc(st.filter)}">` : ""}
        ${options.tools || ""}
        <span class="muted" style="margin-left:auto">${all.length.toLocaleString()} rows</span>
      </div>
      <div class="tbl-wrap"><table class="tbl"><thead><tr>${head}</tr></thead>
      <tbody>${body || `<tr><td colspan="${columns.length}" class="empty">${esc(options.empty || "Nothing to show")}</td></tr>`}</tbody></table></div>
      ${pages > 1 ? `<div class="tbl-foot">
        <button data-page="prev" ${st.page === 0 ? "disabled" : ""}>‹ Prev</button>
        <span>Page ${st.page + 1} of ${pages}</span>
        <button data-page="next" ${st.page >= pages - 1 ? "disabled" : ""}>Next ›</button></div>` : ""}`;

    el.querySelectorAll("th[data-key]").forEach(th => th.addEventListener("click", () => {
      const key = th.dataset.key;
      if (st.sortKey === key) st.sortDir = st.sortDir === "asc" ? "desc" : "asc";
      else { st.sortKey = key; st.sortDir = (columns.find(c => c.key === key) || {}).num ? "desc" : "asc"; }
      st.page = 0;
      draw();
    }));
    el.querySelectorAll("button[data-page]").forEach(b => b.addEventListener("click", () => {
      st.page += b.dataset.page === "next" ? 1 : -1;
      draw();
    }));
    const search = el.querySelector(".tbl-tools input[type=search]");
    if (search) {
      search.addEventListener("input", () => { st.filter = search.value; st.page = 0; draw(); keepFocus(); });
    }
    if (options.onDraw) options.onDraw(el);
  }

  function keepFocus() {
    const input = el.querySelector(".tbl-tools input[type=search]");
    if (input) { input.focus(); input.setSelectionRange(input.value.length, input.value.length); }
  }

  draw();
  return { redraw: draw, setRows(next) { rows = next; st.page = 0; draw(); } };
}
