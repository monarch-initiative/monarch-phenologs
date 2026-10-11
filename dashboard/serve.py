"""
Serves dashboard/site locally. A plain static server with HTTP range support, so the browser side
can read parquet the same way it would from a bucket.

    python dashboard/serve.py [--port 8765] [--dir dashboard/site]
"""
import os
import re
import argparse
import http.server


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    extensions_map = dict(http.server.SimpleHTTPRequestHandler.extensions_map,
                          **{".wasm": "application/wasm", ".mjs": "text/javascript", ".js": "text/javascript",
                             ".parquet": "application/octet-stream", ".json": "application/json"})

    def send_head(self):
        self.range = None
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
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--dir", type=str, default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "site"))
    args = parser.parse_args()

    os.chdir(args.dir)
    print("- Serving {} at http://localhost:{}/".format(args.dir, args.port))
    http.server.ThreadingHTTPServer(("", args.port), RangeHandler).serve_forever()
