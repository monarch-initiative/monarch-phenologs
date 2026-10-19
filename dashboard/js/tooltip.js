// One floating tooltip for every chart
let el = null;

function node() {
  if (!el) {
    el = document.createElement("div");
    el.className = "tooltip";
    el.hidden = true;
    document.body.appendChild(el);
  }
  return el;
}

export function showTip(event, html) {
  const t = node();
  t.innerHTML = html;
  t.hidden = false;
  moveTip(event);
}

export function moveTip(event) {
  const t = node();
  const pad = 14;
  let x = event.clientX + pad, y = event.clientY + pad;
  const r = t.getBoundingClientRect();
  if (x + r.width > window.innerWidth - 8) x = event.clientX - r.width - pad;
  if (y + r.height > window.innerHeight - 8) y = event.clientY - r.height - pad;
  t.style.left = x + "px";
  t.style.top = y + "px";
}

export function hideTip() {
  if (el) el.hidden = true;
}
