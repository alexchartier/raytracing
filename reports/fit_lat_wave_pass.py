"""Fit a shared latitude wave to the 20 O/X ionograms, without density truth.

This is a fast, plane-stratified initialization for the expensive 3-D forward
check. It fits the mode-resolved noses and the latitude variation of group
range. Frequency-specific range offsets absorb the systematic difference
between a plane-stratified approximation and PyLap's magnetoionic rays.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

from analyze_d_inverse_fit import _fitted_profile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
PRIOR = DATA / "d_inverse_iri_peak_summary.json"
BACKGROUND = DATA / "d_inverse_iri_pyiri_background.npz"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
OBSERVED = DATA / "lat_wave_pass_ionograms_recovered"
OUTPUT = DATA / "lat_wave_pass_initial_fit.json"


def plane_group_range(alt: np.ndarray, density: np.ndarray, frequency: float) -> float:
    """Two-way ordinary-wave group path for a piecewise-linear N_e profile."""
    plasma_squared = (0.00898 ** 2) * density
    peak = int(np.argmax(density))
    q = frequency * frequency
    if q >= plasma_squared[peak] or q <= plasma_squared[-1]:
        return np.nan
    index = peak + int(np.flatnonzero(plasma_squared[peak:] <= q)[0])
    turning_alt = alt[index - 1] + (q - plasma_squared[index - 1]) * (
        alt[index] - alt[index - 1]) / (plasma_squared[index] - plasma_squared[index - 1])
    heights = np.r_[turning_alt, alt[index:]]
    fractions = np.sqrt(np.clip(1.0 - np.r_[q, plasma_squared[index:]] / q, 0.0, 1.0))
    return float(4.0 * np.sum(np.diff(heights) / np.maximum(
        fractions[1:] + fractions[:-1], 1e-8)))


def load_observations() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                 list[tuple[int, int, float, float]]]:
    # The manifest is used only for coordinates; its wave field is never read.
    manifest = json.loads(MANIFEST.read_text())
    positions = manifest["profiles"]
    prior = json.loads(PRIOR.read_text())["retrieved"]
    latitudes = np.array([p["latitude_deg"] for p in positions])
    profiles = np.array([_fitted_profile(BACKGROUND, p["latitude_deg"],
                                         p["longitude_deg"], prior)[1]
                         for p in positions])
    altitude = np.load(BACKGROUND, allow_pickle=False)["altitudes_km"]
    noses = np.full((len(positions), 2), np.nan)
    ranges = []
    for row, position in enumerate(positions):
        path = OBSERVED / f"ionogram_{position['index']:02d}.npz"
        with np.load(path, allow_pickle=False) as data:
            records = np.asarray(data["records"], dtype=float)
            frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
        for mode_index, mode in enumerate((1, -1)):
            group = records[records[:, 1] == mode]
            if len(group):
                noses[row, mode_index] = frequencies[int(np.max(group[:, 0]))]
            for frequency_index in np.unique(group[:, 0]).astype(int):
                f = frequencies[frequency_index]
                if 2.4 <= f <= (4.6 if mode == 1 else 5.0):
                    values = group[group[:, 0] == frequency_index, 2]
                    ranges.append((row, mode_index, round(float(f), 1),
                                   float(np.median(values))))
    return latitudes, altitude, profiles, noses, ranges


def fit() -> dict:
    latitudes, altitude, background, noses, ranges = load_observations()
    latitude_center = float(np.mean(latitudes))
    x_km = 6371.0088 * np.deg2rad(latitudes - latitude_center)
    x_normal = x_km / max(abs(x_km))
    peak_altitude = altitude[np.argmax(background, axis=1)]
    envelope = np.exp(-0.5 * ((altitude[None, :] - peak_altitude[:, None]) / 65.0) ** 2)
    range_rows = np.array([r[0] for r in ranges], dtype=int)
    range_modes = np.array([r[1] for r in ranges], dtype=int)
    range_frequencies = np.array([r[2] for r in ranges])
    observed_ranges = np.array([r[3] for r in ranges])
    groups = sorted(set(zip(range_modes, range_frequencies)))
    group_index = np.array([groups.index((m, f))
                            for m, f in zip(range_modes, range_frequencies)])
    group_members = [np.flatnonzero(group_index == i) for i in range(len(groups))]
    # Discard frequencies seen at fewer than seven latitudes; they cannot
    # constrain a smooth wave once a mode/frequency offset is removed.
    valid_groups = [members for members in group_members if len(members) >= 7]
    all_indices = np.concatenate(valid_groups)
    range_rows = range_rows[all_indices]
    range_modes = range_modes[all_indices]
    range_frequencies = range_frequencies[all_indices]
    observed_ranges = observed_ranges[all_indices]
    group_members = []
    for m, f in sorted(set(zip(range_modes, range_frequencies))):
        group_members.append(np.flatnonzero((range_modes == m) & (range_frequencies == f)))

    # Parameters: global peak scale, linear background drift over the pass,
    # wave fraction, wavelength, phase, and O/X cutoff offsets in MHz.
    bounds = [(0.90, 1.12), (-0.12, 0.12), (0.0, 0.35),
              (550.0, 1500.0), (-np.pi, np.pi), (-0.05, 0.20), (0.15, 0.65)]
    def field(p: np.ndarray) -> np.ndarray:
        scale, slope, amplitude, wavelength, phase = p[:5]
        wave = np.cos(2 * np.pi * x_km / wavelength + phase)
        factor = scale * (1.0 + slope * x_normal[:, None]
                          + amplitude * wave[:, None] * envelope)
        return background * factor

    def loss(p: np.ndarray, components: bool = False):
        density = field(p)
        peaks = 0.00898 * np.sqrt(np.max(density, axis=1))
        predicted_noses = peaks[:, None] + p[None, 5:7]
        residual = noses - predicted_noses
        # Missing homing near the cutoff pushes observed noses downward.
        # Preserve information from them at lower weight rather than treating
        # every last accepted frequency as an exact critical frequency.
        weight = np.where(residual < 0, 0.35, 1.0)
        weight[:, 0] *= 0.7
        nose_loss = np.mean(weight * np.minimum(residual ** 2, 0.55 ** 2)) / 0.20 ** 2
        predictions = np.empty(len(range_rows))
        cache = {}
        for i, (row, mode, frequency) in enumerate(zip(range_rows, range_modes, range_frequencies)):
            key = (row, mode, frequency)
            if key not in cache:
                effective_frequency = frequency - (p[6] if mode == 1 else p[5])
                cache[key] = plane_group_range(altitude, density[row], effective_frequency)
            predictions[i] = cache[key]
        residual_ranges = observed_ranges - predictions
        range_terms = []
        for indices in group_members:
            values = residual_ranges[indices]
            values = values[np.isfinite(values)]
            if len(values) >= 7:
                centered = values - np.median(values)
                range_terms.extend(np.minimum((centered / 65.0) ** 2, 4.0))
        range_loss = float(np.mean(range_terms)) if range_terms else 10.0
        penalty = ((p[0] - 1) / .08) ** 2 * .03 + (p[1] / .08) ** 2 * .03
        total = float(nose_loss + 0.75 * range_loss + penalty)
        return (total, nose_loss, range_loss, penalty) if components else total

    solution = differential_evolution(loss, bounds, seed=20260927,
                                      popsize=9, maxiter=28, polish=False,
                                      workers=1, updating="immediate")
    polished = minimize(loss, solution.x, bounds=bounds, method="Powell",
                        options={"maxiter": 80, "xtol": 2e-3, "ftol": 2e-3})
    chosen = polished.x if polished.fun < solution.fun else solution.x
    density = field(chosen)
    modeled_noses = (0.00898 * np.sqrt(np.max(density, axis=1)))[:, None] + chosen[None, 5:7]
    result = {
        "method": "joint O/X noses and centered O/X group-range latitude variation",
        "selection_uses_truth_density": False,
        "background": str(BACKGROUND.relative_to(ROOT)),
        "observed_ionograms": str(OBSERVED.relative_to(ROOT)),
        "global_density_scale": float(chosen[0]),
        "linear_background_fraction_end_to_end_halfspan": float(chosen[1]),
        "wave_amplitude_fraction_at_local_f2_peak": float(chosen[2]),
        "wave_wavelength_km": float(chosen[3]),
        "wave_phase_rad_at_pass_center": float(chosen[4]),
        "wave_vertical_sigma_km_assumed": 65.0,
        "o_nose_offset_mhz": float(chosen[5]),
        "x_nose_offset_mhz": float(chosen[6]),
        "loss": dict(zip(("total", "nose", "range", "prior"), loss(chosen, True))),
        "observed_noses_mhz": noses.tolist(),
        "modeled_noses_mhz": modeled_noses.tolist(),
        "fitted_peak_cm3": np.max(density, axis=1).tolist(),
        "pass_latitudes_deg": latitudes.tolist(),
        "range_samples": len(range_rows),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("observed_noses_mhz", "modeled_noses_mhz",
                                     "fitted_peak_cm3", "pass_latitudes_deg")}, indent=2))
    return result


if __name__ == "__main__":
    fit()
