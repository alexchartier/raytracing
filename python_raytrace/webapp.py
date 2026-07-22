from __future__ import annotations

import argparse
import datetime as dt
import json
import mimetypes
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Iterable
from wsgiref.simple_server import make_server

from .geometry import GeoPoint
from .tracer import PointToPointRayTracer


STATIC_DIR = Path(__file__).resolve().parent / "web_assets"


@dataclass(frozen=True)
class WebTraceConfig:
    alt_min_km: float = 60.0
    alt_max_km: float = 500.0
    alt_step_km: float = 10.0
    lat_step_deg: float = 2.0
    lon_step_deg: float = 2.0
    lat_margin_deg: float = 3.0
    lon_margin_deg: float = 3.0
    nhops: int = 2
    homing_tolerance_m: float = 10000.0
    d_region_model: str = "fpt2018"
    ap_daily: float | None = None
    f107: float | None = None
    refresh_indices: bool = False


def _json_response(start_response, status: str, payload: dict) -> list[bytes]:
    body = json.dumps(payload, indent=2).encode("utf-8")
    headers = [
        ("Content-Type", "application/json; charset=utf-8"),
        ("Content-Length", str(len(body))),
        ("Cache-Control", "no-store"),
    ]
    start_response(status, headers)
    return [body]


def _text_response(start_response, status: str, text: str, content_type: str = "text/plain; charset=utf-8") -> list[bytes]:
    body = text.encode("utf-8")
    start_response(status, [("Content-Type", content_type), ("Content-Length", str(len(body)))])
    return [body]


def _serve_static(start_response, relative_path: str) -> list[bytes]:
    path = (STATIC_DIR / relative_path).resolve()
    static_root = STATIC_DIR.resolve()
    if not path.is_file() or (static_root not in path.parents and path != static_root):
        return _text_response(start_response, "404 Not Found", "Not found")
    body = path.read_bytes()
    content_type, _ = mimetypes.guess_type(str(path))
    start_response(
        "200 OK",
        [
            ("Content-Type", content_type or "application/octet-stream"),
            ("Content-Length", str(len(body))),
        ],
    )
    return [body]


def _read_json_body(environ) -> dict:
    content_length = int(environ.get("CONTENT_LENGTH") or "0")
    raw = environ["wsgi.input"].read(content_length) if content_length > 0 else b"{}"
    return json.loads(raw.decode("utf-8"))


def _parse_iso_datetime(date_text: str, time_text: str) -> dt.datetime:
    parsed_date = dt.date.fromisoformat(date_text)
    parsed_time = dt.time.fromisoformat(time_text)
    return dt.datetime.combine(parsed_date, parsed_time)


def _require_float(payload: dict, key: str) -> float:
    if key not in payload:
        raise ValueError(f"missing field `{key}`")
    return float(payload[key])


def _optional_float(payload: dict, key: str, default: float | None = None) -> float | None:
    value = payload.get(key, default)
    if value in (None, ""):
        return None
    return float(value)


def _trace_payload_to_inputs(payload: dict) -> tuple[dt.datetime, GeoPoint, GeoPoint, float, float | None]:
    when = _parse_iso_datetime(str(payload["date"]), str(payload["time"]))
    tx = GeoPoint(
        _require_float(payload, "tx_lat"),
        _require_float(payload, "tx_lon"),
        float(payload.get("tx_alt_km", 0.0)),
    )
    rx = GeoPoint(
        _require_float(payload, "rx_lat"),
        _require_float(payload, "rx_lon"),
        float(payload.get("rx_alt_km", 0.0)),
    )
    frequency_mhz = _require_float(payload, "frequency_mhz")
    f107 = _optional_float(payload, "f107")
    return when, tx, rx, frequency_mhz, f107


def _trace_result_payload(ray, when: dt.datetime, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float, config: WebTraceConfig, f107: float | None) -> dict:
    absorption_db = None if ray.total_absorption_db is None else float(ray.total_absorption_db)
    return {
        "input": {
            "date_utc": when.date().isoformat(),
            "time_utc": when.time().isoformat(timespec="minutes"),
            "tx": asdict(tx),
            "rx": asdict(rx),
            "frequency_mhz": float(frequency_mhz),
            "f107": None if f107 is None else float(f107),
        },
        "config": asdict(config),
        "result": {
            "reachable": bool(ray.home),
            "error_m": float(ray.error_m),
            "absorption_db": absorption_db,
            "launch_elevation_deg": float(ray.launch_elevation_deg),
            "launch_bearing_deg": float(ray.launch_bearing_deg),
            "group_range_km": None if ray.group_range_to_rx_km is None else float(ray.group_range_to_rx_km),
            "geometric_distance_km": None if ray.geometric_dist_to_rx_km is None else float(ray.geometric_dist_to_rx_km),
            "perigee_km": None if ray.perigee_km is None else float(ray.perigee_km),
            "message": (
                "A converged ray path reached the receiver within the search tolerance."
                if ray.home
                else "No converged ray path reached the receiver within the search tolerance."
            ),
        },
    }


def create_app(
    *,
    tracer_factory: Callable[[], PointToPointRayTracer] | None = None,
    config: WebTraceConfig | None = None,
):
    tracer_factory = tracer_factory or PointToPointRayTracer
    config = config or WebTraceConfig()

    def app(environ, start_response):
        path = environ.get("PATH_INFO", "/")
        method = environ.get("REQUEST_METHOD", "GET").upper()

        if method == "GET" and path == "/":
            return _serve_static(start_response, "index.html")
        if method == "GET" and path.startswith("/static/"):
            return _serve_static(start_response, path.removeprefix("/static/"))
        if method == "GET" and path == "/api/config":
            return _json_response(start_response, "200 OK", {"config": asdict(config)})
        if method == "POST" and path == "/api/trace":
            try:
                payload = _read_json_body(environ)
                when, tx, rx, frequency_mhz, f107 = _trace_payload_to_inputs(payload)
                tracer = tracer_factory()
                ray = tracer.trace_frequencies(
                    when=when,
                    tx=tx,
                    rx=rx,
                    frequencies_mhz=[frequency_mhz],
                    f107=f107,
                    nhops=config.nhops,
                    homing_tolerance_m=config.homing_tolerance_m,
                    alt_min_km=config.alt_min_km,
                    alt_max_km=config.alt_max_km,
                    alt_step_km=config.alt_step_km,
                    lat_step_deg=config.lat_step_deg,
                    lon_step_deg=config.lon_step_deg,
                    lat_margin_deg=config.lat_margin_deg,
                    lon_margin_deg=config.lon_margin_deg,
                    d_region_model=config.d_region_model,
                    ap_daily=config.ap_daily,
                    refresh_indices=config.refresh_indices,
                )[0]
            except KeyError as exc:
                return _json_response(start_response, "400 Bad Request", {"error": f"missing field `{exc.args[0]}`"})
            except ValueError as exc:
                return _json_response(start_response, "400 Bad Request", {"error": str(exc)})
            except Exception as exc:
                return _json_response(start_response, "500 Internal Server Error", {"error": f"raytrace failed: {exc}"})

            return _json_response(
                start_response,
                "200 OK",
                _trace_result_payload(ray, when, tx, rx, frequency_mhz, config, f107),
            )

        return _text_response(start_response, "404 Not Found", "Not found")

    return app


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Prototype HF link visualizer web app.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(list(argv) if argv is not None else None)

    app = create_app()
    with make_server(args.host, args.port, app) as server:
        print(f"Serving python_raytrace visualizer on http://{args.host}:{args.port}")
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
