"""Freeze a public NeQuick-G model profile, separate from all retrieval inputs.

The external source is tpl2go/NequickG at commit 1d1783412ed9811540b56d409a1dcf27d2413740.
Its Python 2 files need a mechanical 2to3 conversion in an ignored cache directory.
This script imports only that external model and NumPy; it does not import any
retrieval or ray-tracing modules.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/data/nequick_independent_case"
COMMIT = "1d1783412ed9811540b56d409a1dcf27d2413740"
LATITUDE = -53.67859907285381
LONGITUDE = 7.7400542554165135
MONTH = 1
UTC_HOUR = 12.0
AZ = 64.0
LAT_GRID = np.arange(-80.0, -28.0 + 0.1, 2.0)
LON_GRID = np.arange(-24.0, 40.0 + 0.1, 4.0)
ALT_GRID = np.arange(60.0, 900.0 + 0.1, 20.0)


def main(source: Path, out: Path = OUT, latitude: float = LATITUDE,
         longitude: float = LONGITUDE, month: int = MONTH,
         utc_hour: float = UTC_HOUR, az: float = AZ) -> None:
    source = source.resolve()
    if not (source / "NequickG.py").is_file():
        raise FileNotFoundError(source / "NequickG.py")
    sys.path.insert(0, str(source))
    from NequickG import GalileoBroadcast, NEQTime, Position
    from NequickG_global import NequickG_global

    model = NequickG_global(NEQTime(month, utc_hour),
                           GalileoBroadcast(az, 0.0, 0.0))
    local, parameters = model.get_Nequick_local(Position(latitude, longitude))
    high_alt = np.arange(60.0, 900.0 + 0.1, 1.0)
    high_density = np.asarray(local.electrondensity(high_alt), dtype=float) / 1e6
    sampled_density = np.asarray(local.electrondensity(ALT_GRID), dtype=float) / 1e6
    if (not np.all(np.isfinite(high_density)) or np.any(high_density <= 0)
            or not np.all(np.isfinite(sampled_density)) or np.any(sampled_density <= 0)):
        raise ValueError("NeQuick-G produced nonpositive or nonfinite density")
    grid = np.broadcast_to(sampled_density,
                           (len(LAT_GRID), len(LON_GRID), len(ALT_GRID))).copy()
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out / "truth_density.npz", latitudes_deg=LAT_GRID,
        longitudes_deg=LON_GRID, altitudes_km=ALT_GRID,
        electron_density_cm3=grid,
        model=np.array("NeQuick-G; tpl2go implementation"),
        time_utc=np.array(f"month {month:02d}, {utc_hour:g} UT"))
    np.savez_compressed(out / "truth_profile_1km.npz",
                        altitudes_km=high_alt, electron_density_cm3=high_density)
    metadata = {
        "model": "NeQuick-G, Galileo ionospheric correction model",
        "public_source": "https://github.com/tpl2go/NequickG",
        "source_commit": COMMIT,
        "conversion": "Python 2 to Python 3 via lib2to3; no model-equation edits",
        "converted_nequickg_py_sha256": hashlib.sha256(
            (source / "NequickG.py").read_bytes()).hexdigest(),
        "latitude_deg": latitude,
        "longitude_deg": longitude,
        "month": month,
        "utc_hour": utc_hour,
        "galileo_broadcast_ai0": az,
        "galileo_broadcast_ai1": 0.0,
        "galileo_broadcast_ai2": 0.0,
        "model_fof2_mhz": float(parameters.foF2),
        "model_hmf2_km": float(parameters.hmF2),
        "high_resolution_peak_altitude_km": float(high_alt[np.argmax(high_density)]),
        "high_resolution_peak_density_cm3": float(np.max(high_density)),
        "horizontal_field": "one NeQuick-G profile replicated across ray-tracing grid",
    }
    (out / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="Pinned external NeQuick-G checkout after 2to3 conversion")
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--latitude", type=float, default=LATITUDE)
    parser.add_argument("--longitude", type=float, default=LONGITUDE)
    parser.add_argument("--month", type=int, default=MONTH)
    parser.add_argument("--utc-hour", type=float, default=UTC_HOUR)
    parser.add_argument("--az", type=float, default=AZ)
    args = parser.parse_args()
    main(args.source, args.out, args.latitude, args.longitude,
         args.month, args.utc_hour, args.az)
