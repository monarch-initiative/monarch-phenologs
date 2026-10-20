"""
Serves the dashboard locally, with a run's releases at /releases/ (where the site looks for them) and a
releases.json listing them. A plain static server with HTTP range support, so the browser reads parquet the
way it would from the bucket.

    python dashboard/serve.py -p path/to/project_dir [--port 8765]
"""
import io
import os
import re
import glob
import json
import argparse
import http.server


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    releases = None
    extensions_map = dict(http.server.SimpleHTTPRequestHandler.extensions_map,
                          **{".wasm": "application/wasm", ".mjs": "text/javascript", ".js": "text/javascript",
                             ".parquet": "application/octet-stream", ".json": "application/json"})

    def translate_path(self, path):
        """ /releases/... is read from the run's release folder, everything else from the dashboard """
        site = super().translate_path(path)
        mount = os.path.join(self.directory, "releases")
        if site == mount or site.startswith(mount + os.sep):
            return self.releases + site[len(mount):]
        return site

    def send_head(self):
        self.range = None
        if self.path.split("?")[0] == "/releases/releases.json":
            return self.release_list()
        header = self.headers.get("Range")
        m = re.match(r"bytes=(\d*)-(\d*)$", header or "")
        path = self.translate_path(self.path)
        if m and os.path.isfile(path):
            size = os.path.getsize(path)
            start = int(m.group(1)) if m.group(1) else max(0, size - int(m.group(2)))
            end = min(int(m.group(2)), size - 1) if m.group(1) and m.group(2) else size - 1
            if start <= end:
                self.range = (start, end)
                fh = open(path, "rb")
                self.send_response(206)
                self.send_header("Content-Type", self.guess_type(path))
                self.send_header("Content-Range", "bytes {}-{}/{}".format(start, end, size))
                self.send_header("Content-Length", str(end - start + 1))
                self.send_header("Accept-Ranges", "bytes")
                self.end_headers()
                fh.seek(start)
                return fh
        return super().send_head()

    def release_list(self):
        """ The run's releases, listed the way the published releases.json lists them """
        entries = []
        for path in glob.glob(os.path.join(self.releases, "*", "*", "overview.json")):
            with open(path) as fh:
                entries.append(dict(json.load(fh)["release"], path=os.path.relpath(os.path.dirname(path), self.releases) + "/"))
        body = json.dumps(sorted(entries, key=lambda r: r["graph_version"], reverse=True)).encode()
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        return io.BytesIO(body)

    def copyfile(self, source, outputfile):
        if self.range is None:
            return super().copyfile(source, outputfile)
        remaining = self.range[1] - self.range[0] + 1
        while remaining > 0:
            chunk = source.read(min(1 << 20, remaining))
            if not chunk:
                break
            outputfile.write(chunk)
            remaining -= len(chunk)

    def end_headers(self):
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, fmt, *args):
        pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Serve the phenologs dashboard locally")
    parser.add_argument("-p", "--project_dir", help="Project directory of a run built with python/07_build_release.py", required=True, type=str)
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()

    RangeHandler.releases = os.path.join(os.path.abspath(args.project_dir), "release")
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print("- Serving the dashboard at http://localhost:{}/ with the releases in {}".format(args.port, RangeHandler.releases))
    http.server.ThreadingHTTPServer.request_queue_size = 64     # the browser asks for every module at once; the default 5 drops some
    http.server.ThreadingHTTPServer(("", args.port), RangeHandler).serve_forever()
