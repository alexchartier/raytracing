"""Build symmetric nearby pass candidates before Doppler-based selection.

The candidate perturbations are fixed from the ionogram-only solution. This
script does not open the withheld density truth.
"""

from __future__ import annotations

import json
from pathlib import Path

from build_lat_wave_retrieval_grid import build

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
BASE = DATA / "lat_wave_pass_initial_fit.json"
OUTPUT = DATA / "lat_wave_doppler_candidates"
PARAMETERS = {
    "scale_low": ("global_density_scale", -0.02),
    "scale_high": ("global_density_scale", 0.02),
    "amplitude_low": ("wave_amplitude_fraction_at_local_f2_peak", -0.03),
    "amplitude_high": ("wave_amplitude_fraction_at_local_f2_peak", 0.03),
}


def main() -> None:
    original = json.loads(BASE.read_text())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    candidates = []
    for name, (parameter, offset) in PARAMETERS.items():
        fit = dict(original)
        fit[parameter] = original[parameter] + offset
        fit["doppler_candidate_change"] = {"parameter": parameter, "offset": offset}
        fit_path = OUTPUT / f"{name}.json"
        density_path = OUTPUT / f"{name}_density.npz"
        grid_path = OUTPUT / f"{name}_grid.nc"
        fit_path.write_text(json.dumps(fit, indent=2) + "\n")
        build(fit_path, density_path, grid_path)
        candidates.append({"name": name, "fit": str(fit_path.relative_to(ROOT)),
                           "density": str(density_path.relative_to(ROOT)),
                           "grid": str(grid_path.relative_to(ROOT)),
                           "parameter": parameter, "offset": offset})
    (OUTPUT / "candidates.json").write_text(json.dumps(candidates, indent=2) + "\n")


if __name__ == "__main__":
    main()
