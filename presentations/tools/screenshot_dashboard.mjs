// Screenshot a page of the phenologs dashboard for the "resource" slide.
//   node tools/screenshot_dashboard.mjs [url] [out.png]
// The site loads its tables with DuckDB-WASM, so a plain `chrome --screenshot` fires too early.
// This drives headless Chrome over the DevTools protocol and waits for the page's own text.
// Needs Google Chrome and Node 22+ (built-in WebSocket). Uses a throwaway profile.
import { spawn } from "node:child_process";
import { mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const url = process.argv[2] || "https://monarch-initiative.github.io/monarch-phenologs/?network=disease#/disease/MONDO:0015229";
const out = resolve(process.argv[3] || join(here, "..", "assets", "dashboard_bbs.png"));
const WAIT_FOR = process.env.WAIT_FOR || "CANDIDATE GENES";      // text that means the page has its data
const [W, H] = (process.env.VIEWPORT || "1500x900").split("x").map(Number);
const CHROME = process.env.CHROME || "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const PORT = 9300 + Math.floor(Math.random() * 500);

const profile = mkdtempSync(join(tmpdir(), "phenologs-shot-"));
const chrome = spawn(CHROME, ["--headless=new", `--remote-debugging-port=${PORT}`, `--user-data-dir=${profile}`,
  "--no-first-run", "--hide-scrollbars", "--disable-gpu", "about:blank"], { stdio: "ignore" });
const sleep = ms => new Promise(r => setTimeout(r, ms));

async function devtools() {
  for (let i = 0; i < 50; i++) {
    try { return (await (await fetch(`http://127.0.0.1:${PORT}/json/version`)).json()).webSocketDebuggerUrl; }
    catch { await sleep(200); }
  }
  throw new Error("Chrome did not start");
}

try {
  const ws = new WebSocket(await devtools());
  await new Promise((ok, bad) => { ws.onopen = ok; ws.onerror = bad; });
  let id = 0;
  const waiting = new Map();
  ws.onmessage = e => { const m = JSON.parse(e.data); if (m.id && waiting.has(m.id)) { waiting.get(m.id)(m); waiting.delete(m.id); } };
  const send = (method, params = {}, sessionId) => new Promise((ok, bad) => {
    const n = ++id;
    waiting.set(n, m => m.error ? bad(new Error(`${method}: ${m.error.message}`)) : ok(m.result));
    ws.send(JSON.stringify({ id: n, method, params, sessionId }));
  });

  const { targetId } = await send("Target.createTarget", { url: "about:blank" });
  const { sessionId } = await send("Target.attachToTarget", { targetId, flatten: true });
  await send("Emulation.setDeviceMetricsOverride", { width: W, height: H, deviceScaleFactor: 2, mobile: false }, sessionId);
  await send("Page.enable", {}, sessionId);
  await send("Page.navigate", { url }, sessionId);
  let ready = false;
  for (let i = 0; i < 120 && !ready; i++) {
    await sleep(500);
    const r = await send("Runtime.evaluate", { expression: `document.body.innerText.includes(${JSON.stringify(WAIT_FOR)})`, returnByValue: true }, sessionId);
    ready = r.result.value === true;
  }
  if (!ready) throw new Error(`"${WAIT_FOR}" never appeared at ${url}`);
  await sleep(1500);                                              // let the charts finish drawing
  const { data } = await send("Page.captureScreenshot", { format: "png" }, sessionId);
  writeFileSync(out, Buffer.from(data, "base64"));
  console.log(out);
  ws.close();
} finally {
  chrome.kill();
  await sleep(300);
  rmSync(profile, { recursive: true, force: true });
}
