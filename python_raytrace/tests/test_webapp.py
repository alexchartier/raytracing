import io
import json
import unittest

from python_raytrace.geometry import GeoPoint
from python_raytrace.webapp import WebTraceConfig, create_app


class _FakeRay:
    home = True
    error_m = 1234.5
    total_absorption_db = 9.75
    group_range_to_rx_km = 1800.0
    geometric_dist_to_rx_km = 1700.0
    perigee_km = 0.0

    @property
    def launch_elevation_deg(self):
        return 22.0

    @property
    def launch_bearing_deg(self):
        return 145.0


class _FakeTracer:
    def __init__(self):
        self.calls = []

    def trace_frequencies(self, **kwargs):
        self.calls.append(kwargs)
        return [_FakeRay()]


def _invoke_app(app, *, method: str, path: str, body: bytes = b"", content_type: str = "application/json"):
    status_headers = {}

    def start_response(status, headers):
        status_headers["status"] = status
        status_headers["headers"] = dict(headers)

    environ = {
        "REQUEST_METHOD": method,
        "PATH_INFO": path,
        "CONTENT_LENGTH": str(len(body)),
        "CONTENT_TYPE": content_type,
        "wsgi.input": io.BytesIO(body),
    }
    response_body = b"".join(app(environ, start_response))
    return status_headers["status"], status_headers["headers"], response_body


class WebAppTests(unittest.TestCase):
    def test_index_serves_html(self) -> None:
        app = create_app()
        status, headers, body = _invoke_app(app, method="GET", path="/")
        self.assertEqual(status, "200 OK")
        self.assertIn("text/html", headers["Content-Type"])
        self.assertIn(b"HF Path Visualizer", body)

    def test_trace_endpoint_returns_json(self) -> None:
        tracer = _FakeTracer()
        app = create_app(tracer_factory=lambda: tracer, config=WebTraceConfig())
        payload = {
            "date": "2020-01-15",
            "time": "12:00",
            "tx_lat": -77.8,
            "tx_lon": 166.4,
            "tx_alt_km": 0.0,
            "rx_lat": -89.9,
            "rx_lon": 166.4,
            "rx_alt_km": 1.0,
            "frequency_mhz": 4.1,
        }
        status, headers, body = _invoke_app(
            app,
            method="POST",
            path="/api/trace",
            body=json.dumps(payload).encode("utf-8"),
        )
        self.assertEqual(status, "200 OK")
        self.assertIn("application/json", headers["Content-Type"])
        parsed = json.loads(body.decode("utf-8"))
        self.assertTrue(parsed["result"]["reachable"])
        self.assertEqual(parsed["result"]["absorption_db"], 9.75)
        self.assertEqual(tracer.calls[0]["alt_max_km"], 500.0)
        self.assertEqual(tracer.calls[0]["nhops"], 2)
        self.assertIsNone(tracer.calls[0]["f107"])

    def test_trace_endpoint_validates_missing_fields(self) -> None:
        app = create_app(tracer_factory=_FakeTracer)
        status, _, body = _invoke_app(
            app,
            method="POST",
            path="/api/trace",
            body=json.dumps({"date": "2020-01-15"}).encode("utf-8"),
        )
        self.assertEqual(status, "400 Bad Request")
        parsed = json.loads(body.decode("utf-8"))
        self.assertIn("missing field", parsed["error"])


if __name__ == "__main__":
    unittest.main()
