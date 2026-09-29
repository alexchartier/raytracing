"""Freeze an 800 km in situ density observation for the April NeQuick case.

Only this preparation step reads the independent profile. The fitting code
receives the scalar observation and the O/X ionogram, not the profile.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "reports/data/nequick_prospective_case"
DESTINATION = ROOT / "reports/data/nequick_local_density_case"


def main() -> None:
    DESTINATION.mkdir(parents=True, exist_ok=True)
    for name in ("truth_ionogram_recovered.npz", "truth_profile_1km.npz",
                 "truth_density.npz", "provenance.json"):
        shutil.copyfile(SOURCE / name, DESTINATION / name)
    with np.load(SOURCE / "truth_profile_1km.npz", allow_pickle=False) as data:
        altitude = np.asarray(data["altitudes_km"], dtype=float)
        density = np.asarray(data["electron_density_cm3"], dtype=float)
    spacecraft_altitude = 800.0
    observation = {
        "spacecraft_altitude_km": spacecraft_altitude,
        "electron_density_cm3": float(np.interp(spacecraft_altitude, altitude, density)),
        "measurement_type": "assumed exact local in situ electron density",
        "source": "April independent NeQuick-G truth profile at spacecraft altitude",
        "use_during_fit": "one scalar only; no other truth density or peak parameter",
    }
    (DESTINATION / "local_density_observation.json").write_text(
        json.dumps(observation, indent=2) + "\n")
    print(json.dumps(observation, indent=2))


if __name__ == "__main__":
    main()
