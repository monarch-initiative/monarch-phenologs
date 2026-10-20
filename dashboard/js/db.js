// DuckDB-WASM over a release's parquet files. Tables are fetched the first time a page needs them
// and stay registered for the rest of the session. All SQL the site runs goes through query().
import * as duckdb from "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm";

// Where the releases live: the Monarch data bucket, or releases/ beside the site when it is served locally
// (serve.py). ?data=<url> points it anywhere else
const params = new URLSearchParams(location.search);
const LOCAL = ["localhost", "127.0.0.1"].includes(location.hostname);
export const BUCKET = "https://data.monarchinitiative.org/phenologs/alpha/";
export const ROOT = (params.get("data") || (LOCAL ? "releases" : BUCKET)).replace(/\/?$/, "/");

// The release layout this site reads (SCHEMA in python/phenologs_release_utils.py)
const SCHEMA = 2;

// The two tables meant for download carry their network in the file name (FILE_NAMES in python/phenologs_release_utils.py)
const FILE_NAMES = { phenologs: "{prd}-phenologs", candidates: "gene-{prd}-candidates" };

let releases = [];              // releases.json: every published release, newest first
let open = null;                // the open release's entry in it
const overviews = new Map();    // path -> Promise of overview.json
let db = null;
let conn = null;
const registered = new Map();   // table -> Promise

// Nothing to open: no release published yet, or none for the KG version or network asked for. The page says so
// instead of showing an error
export class NoRelease extends Error {}

// Opens ?release=<graph_version> (the newest by default) of ?network=disease|phenotype (disease by default)
// from releases.json. Returns its overview.json. The list is asked for with a query no cache has seen, so a
// release shows up as soon as it is published (the bucket otherwise serves copies up to an hour old)
export async function openRelease() {
  const res = await fetch(ROOT + "releases.json?nocache=" + Date.now());
  if (res.status === 404 || res.status === 403) throw new NoRelease("No releases have been published yet.");
  if (!res.ok) throw new Error(`Could not load ${ROOT}releases.json (${res.status})`);
  releases = await res.json();
  if (!releases.length) throw new NoRelease("No releases have been published yet.");

  const want = params.get("release"), network = params.get("network") || "disease";
  open = releases.find(x => x.prediction_network === network && (!want || x.graph_version === want));
  if (!open) throw new NoRelease(`No ${network} release${want ? " for Monarch KG " + want : ""}.`);
  const ov = await overviewOf(open);
  if (ov.release.schema !== SCHEMA) {
    throw new Error(`This release has layout ${ov.release.schema} and this site reads layout ${SCHEMA}. Rebuild the release with the current code`);
  }
  return ov;
}

export function releaseList() {
  return releases;
}

export function openEntry() {
  return open;
}

// A release's overview.json, fetched once
export function overviewOf(entry) {
  if (!overviews.has(entry.path)) overviews.set(entry.path, getJSON(ROOT + entry.path + "overview.json"));
  return overviews.get(entry.path);
}

// The address of a file of a release (the open one by default)
export function fileURL(name, entry = open) {
  return new URL(ROOT + entry.path + name, location.href).href;
}

// The file a table is in
export function fileOf(table, entry = open) {
  return (FILE_NAMES[table] || table).replace("{prd}", entry.prediction_network) + ".parquet";
}

export async function getJSON(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Could not load ${url} (${res.status})`);
  return res.json();
}

async function init() {
  if (conn) return conn;
  const bundle = await duckdb.selectBundle(duckdb.getJsDelivrBundles());
  const workerUrl = URL.createObjectURL(new Blob([`importScripts("${bundle.mainWorker}");`], { type: "text/javascript" }));
  const worker = new Worker(workerUrl);
  db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), worker);
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
  URL.revokeObjectURL(workerUrl);
  conn = await db.connect();
  return conn;
}

// Make the release's file of a table queryable as that table. Fetched whole once
export function ensure(...tables) {
  return Promise.all(tables.map(t => {
    if (!registered.has(t)) registered.set(t, register(t));
    return registered.get(t);
  }));
}

async function register(table) {
  await init();
  const res = await fetch(ROOT + open.path + fileOf(table));
  if (!res.ok) throw new Error(`Could not load ${fileOf(table)} (${res.status})`);
  await db.registerFileBuffer(table + ".parquet", new Uint8Array(await res.arrayBuffer()));
  await conn.query(`CREATE VIEW ${table} AS SELECT * FROM parquet_scan('${table}.parquet')`);
}

// Rows as plain objects; 64 bit integers come back as BigInt from Arrow, so they are made numbers here
export async function query(sql) {
  await init();
  const result = await conn.query(sql);
  return result.toArray().map(row => {
    const o = row.toJSON();
    for (const k in o) if (typeof o[k] === "bigint") o[k] = Number(o[k]);
    return o;
  });
}

// A string literal for SQL
export function lit(value) {
  return "'" + String(value).replace(/'/g, "''") + "'";
}
