from __future__ import annotations

import argparse
import datetime as dt
import json

from .geometry import GeoPoint
from .tracer import PointToPointRayTracer


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Point-to-point HF raytracing using PyIRI, IRI2020 D-region profiles, PyMSIS, and pylap.")
    parser.add_argument("--time", required=True, help="UTC time in ISO-8601 format, for example 2020-01-15T12:00:00")
    parser.add_argument("--tx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"), required=True)
    parser.add_argument("--rx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"), required=True)
    parser.add_argument("--freq", dest="freqs", action="append", type=float, required=True,
                        help="Frequency in MHz. Repeat this flag for a sweep.")
    parser.add_argument("--f107", type=float, default=None, help="Override daily F10.7. If omitted, pymsis indices are used.")
    parser.add_argument("--ap-daily", type=float, default=None, help="Override daily Ap. If omitted, pymsis indices are used.")
    parser.add_argument("--refresh-indices", action="store_true", help="Refresh the cached pymsis F10.7/Ap file before running.")
    parser.add_argument("--ox-mode", type=int, default=0, choices=(-1, 0, 1))
    parser.add_argument("--nhops", type=int, default=2)
    parser.add_argument("--homing-tolerance-m", type=float, default=100.0)
    parser.add_argument("--alt-min-km", type=float, default=60.0)
    parser.add_argument("--alt-max-km", type=float, default=500.0)
    parser.add_argument("--alt-step-km", type=float, default=2.0)
    parser.add_argument("--lat-step-deg", type=float, default=1.0)
    parser.add_argument("--lon-step-deg", type=float, default=1.0)
    parser.add_argument("--lat-margin-deg", type=float, default=5.0)
    parser.add_argument("--lon-margin-deg", type=float, default=5.0)
    parser.add_argument("--d-region-model", choices=("none", "iri1990", "fpt2018"), default="fpt2018")
    parser.add_argument("--blend-bottom-km", type=float, default=120.0)
    parser.add_argument("--blend-top-km", type=float, default=140.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    tracer = PointToPointRayTracer()
    results = tracer.trace_frequencies(
        when=_parse_time(args.time),
        tx=GeoPoint(*args.tx),
        rx=GeoPoint(*args.rx),
        frequencies_mhz=args.freqs,
        f107=args.f107,
        ap_daily=args.ap_daily,
        ox_mode=args.ox_mode,
        nhops=args.nhops,
        homing_tolerance_m=args.homing_tolerance_m,
        alt_min_km=args.alt_min_km,
        alt_max_km=args.alt_max_km,
        alt_step_km=args.alt_step_km,
        lat_step_deg=args.lat_step_deg,
        lon_step_deg=args.lon_step_deg,
        lat_margin_deg=args.lat_margin_deg,
        lon_margin_deg=args.lon_margin_deg,
        d_region_model=args.d_region_model,
        blend_bottom_km=args.blend_bottom_km,
        blend_top_km=args.blend_top_km,
        refresh_indices=args.refresh_indices,
    )
    print(json.dumps([result.to_dict() for result in results], indent=2))


if __name__ == "__main__":
    main()
