"""Local-only review app; Julia owns all mathematical data and all geometry."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit, unquote

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT / "web" / "heegaard"
GENERATED = ROOT / "artifacts" / "heegaard" / "app"
SURFACE_GENERATED = ROOT / "artifacts" / "heegaard" / "full_curves"
LOCK = threading.Lock()

class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(WEB), **kwargs)

    def translate_path(self, path):
        path = unquote(urlsplit(path).path)
        if path.startswith("/surface-generated/"):
            candidate = (SURFACE_GENERATED / path[len("/surface-generated/"):]).resolve()
            if not candidate.is_relative_to(SURFACE_GENERATED.resolve()):
                return str(WEB / "not-found")
            return str(candidate)
        if path.startswith("/surface_") and path.endswith(".html"):
            return str(WEB / "surface.html")
        if path.startswith("/generated/"):
            relative = path[len("/generated/"):]
            candidate = (GENERATED / relative).resolve()
            if not candidate.is_relative_to(GENERATED.resolve()):
                return str(WEB / "not-found")
            return str(candidate)
        if path.startswith("/d_") and path.endswith(".html"):
            return str(WEB / "index.html")
        return super().translate_path(path)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        super().end_headers()

    def reply(self, status, data):
        encoded = json.dumps(data, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)

    def do_POST(self):
        if self.path not in ("/api/render", "/api/convert"):
            return self.reply(404, {"error": "Unknown endpoint"})
        origin = self.headers.get("Origin")
        allowed = {f"http://localhost:{self.server.server_port}", f"http://127.0.0.1:{self.server.server_port}"}
        if origin and origin not in allowed:
            return self.reply(403, {"error": "Local review origin required"})
        try:
            size = int(self.headers.get("Content-Length", "0"))
            if not 0 < size <= 5_000_000:
                return self.reply(413, {"error": "Request must be under 5 MB"})
            data = json.loads(self.rfile.read(size))
            if not isinstance(data, dict) or not isinstance(data.get("diagram"), dict):
                return self.reply(400, {"error": "A .hgd.json diagram is required"})
        except (ValueError, UnicodeError):
            return self.reply(400, {"error": "Invalid JSON"})
        if not LOCK.acquire(blocking=False):
            return self.reply(409, {"error": "Another diagram is rendering. Please wait."})
        try:
            with tempfile.TemporaryDirectory(prefix="heegaard-render-") as directory:
                request, response = Path(directory)/"request.json", Path(directory)/"response.json"
                request.write_text(json.dumps(data), encoding="utf-8")
                script = "convert_request.jl" if self.path == "/api/convert" else "render_request.jl"
                process = subprocess.run([self.server.julia, "--startup-file=no", "--compiled-modules=existing", "--pkgimages=no", "-O0",
                    f"--project={ROOT}", str(ROOT/"scripts/heegaard"/script), str(request), str(response)],
                    cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=900,
                    creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
                if process.returncode != 0:
                    return self.reply(422, {"error": process.stderr[-3500:] or "Julia rendering failed"})
                return self.reply(200, json.loads(response.read_text(encoding="utf-8")))
        except subprocess.TimeoutExpired:
            return self.reply(504, {"error": "Julia processing reached the 15-minute limit; no result was accepted."})
        except OSError as error:
            return self.reply(500, {"error": str(error)})
        finally:
            LOCK.release()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8731)
    parser.add_argument("--julia", default=shutil.which("julia"))
    args = parser.parse_args()
    if not args.julia:
        parser.error("Julia not found; use --julia PATH")
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    server.julia = args.julia
    print(f"Heegaard Diagram: http://localhost:{args.port}", flush=True)
    server.serve_forever()
