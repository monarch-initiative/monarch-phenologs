// DuckDB-WASM over the release's parquet files. Tables are fetched the first time a page needs them
// and stay registered for the rest of the session. All SQL the site runs goes through query().
import * as duckdb from "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm";

const DATA = "data/";
let db = null;
let conn = null;
const registered = new Map();   // table -> Promise

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

// Make data/<table>.parquet queryable as <table>. Fetched whole once; the biggest table is ~15 MB.
export function ensure(...tables) {
  return Promise.all(tables.map(t => {
    if (!registered.has(t)) registered.set(t, register(t));
    return registered.get(t);
  }));
}

async function register(table) {
  await init();
  const res = await fetch(DATA + table + ".parquet");
  if (!res.ok) throw new Error(`Could not load ${table}.parquet (${res.status})`);
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
