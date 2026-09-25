from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.ndimage import gaussian_filter

PACKAGE_ROOT = Path(__file__).resolve().parent.parent
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

from python_raytrace.multisat_topside_inverse_demo import (
    TopsideInverseConfig,
    IonosphereFitParams,
    _apply_fit_params_to_grid,
    _fit_parameter_bounds,
    _line_of_sight_angles,
    _render_frequency_plot_slice,
    _subset_problem_cases,
    build_inverse_problem,
    make_frequency_axis_mhz,
)


DEFAULT_CACHE = PACKAGE_ROOT / ".cache/multisat_topside_inverse/global_background_eeab2b44d517287b.nc"
DEFAULT_CASE = "plane01_sv121_oblique"
DEFAULT_TRUTH = IonosphereFitParams(
    density_scale=0.90,
    hmf2_shift_km=28.0,
    wave_amplitude_fraction=0.22,
    wave_phase_rad=2.35,
    wave_bearing_deg=-42.0,
)


def _wrap_bearing_deg(value: float) -> float:
    return ((float(value) + 180.0) % 360.0) - 180.0


def _clip_params(params: IonosphereFitParams) -> IonosphereFitParams:
    bounds = _fit_parameter_bounds()
    return IonosphereFitParams(
        density_scale=float(np.clip(params.density_scale, *bounds["density_scale"])),
        hmf2_shift_km=float(np.clip(params.hmf2_shift_km, *bounds["hmf2_shift_km"])),
        wave_amplitude_fraction=float(np.clip(params.wave_amplitude_fraction, *bounds["wave_amplitude_fraction"])),
        wave_phase_rad=float(np.clip(params.wave_phase_rad, *bounds["wave_phase_rad"])),
        wave_bearing_deg=_wrap_bearing_deg(np.clip(params.wave_bearing_deg, *bounds["wave_bearing_deg"])),
    )


def _finalize_raw(
    raw: dict[str, np.ndarray],
    config: TopsideInverseConfig,
) -> dict[str, np.ndarray]:
    sigma = (config.blur_sigma_frequency_bins, config.blur_sigma_range_bins)
    o_image = gaussian_filter(np.asarray(raw["o_image"], dtype=float), sigma=sigma, mode="nearest")
    x_image = gaussian_filter(np.asarray(raw["x_image"], dtype=float), sigma=sigma, mode="nearest")
    doppler_num = gaussian_filter(np.asarray(raw["doppler_num"], dtype=float), sigma=sigma, mode="nearest")
    doppler_den = gaussian_filter(np.asarray(raw["doppler_den"], dtype=float), sigma=sigma, mode="nearest")
    total_image = o_image + x_image
    total_scale = max(float(np.max(total_image)), 1e-9)
    o_norm = o_image / total_scale
    x_norm = x_image / total_scale
    total_norm = total_image / total_scale
    ox_split = (o_image - x_image) / np.maximum(total_image, 1e-6)
    doppler_image_hz = np.divide(doppler_num, doppler_den, out=np.zeros_like(doppler_num), where=doppler_den > 1e-9)
    return {
        "o_image": o_norm,
        "x_image": x_norm,
        "total_image": total_norm,
        "ox_split_image": ox_split,
        "doppler_image_hz": doppler_image_hz,
        "doppler_weight": np.clip(total_norm, 0.0, 1.0),
        "range_edges_km": np.asarray(raw["range_edges_km"], dtype=float),
        "range_centers_km": np.asarray(raw["range_centers_km"], dtype=float),
    }


def _single_case_cost(observed: dict[str, np.ndarray], predicted: dict[str, np.ndarray]) -> float:
    doppler_scale_hz = 50.0
    image_cost = float(np.mean((predicted["total_image"] - observed["total_image"]) ** 2))
    split_weight = np.maximum(observed["total_image"], 0.0)
    if float(np.sum(split_weight)) > 0.0:
        split_cost = float(
            np.sum(((predicted["ox_split_image"] - observed["ox_split_image"]) ** 2) * split_weight)
            / np.sum(split_weight)
        )
    else:
        split_cost = 0.0
    doppler_weight = np.maximum(observed["doppler_weight"], 0.0)
    if float(np.sum(doppler_weight)) > 0.0:
        doppler_cost = float(
            np.sum(
                (((predicted["doppler_image_hz"] - observed["doppler_image_hz"]) / doppler_scale_hz) ** 2)
                * doppler_weight
            )
            / np.sum(doppler_weight)
        )
    else:
        doppler_cost = 0.0
    return image_cost + 0.7 * doppler_cost + 0.5 * split_cost


def _frequency_subset(finalized: dict[str, np.ndarray], indices: np.ndarray) -> dict[str, np.ndarray]:
    subset = {key: value for key, value in finalized.items() if key not in {"range_edges_km", "range_centers_km"}}
    return {
        **{key: np.asarray(value, dtype=float)[indices] for key, value in subset.items()},
        "range_edges_km": np.asarray(finalized["range_edges_km"], dtype=float),
        "range_centers_km": np.asarray(finalized["range_centers_km"], dtype=float),
    }


def _candidate_key(params: IonosphereFitParams) -> tuple[float, float, float, float, float]:
    clipped = _clip_params(params)
    return (
        round(clipped.density_scale, 6),
        round(clipped.hmf2_shift_km, 6),
        round(clipped.wave_amplitude_fraction, 6),
        round(clipped.wave_phase_rad, 6),
        round(clipped.wave_bearing_deg, 6),
    )


def _candidate_record(cost: float, params: IonosphereFitParams) -> dict[str, object]:
    return {
        "cost": float(cost),
        "params": vars(_clip_params(params)),
    }


def _child_trace_block(args: argparse.Namespace) -> None:
    frequencies_mhz = tuple(float(value) for value in json.loads(args.child_frequencies_json))
    params = IonosphereFitParams(**json.loads(args.child_params_json))
    config = replace(
        TopsideInverseConfig(),
        planes=(int(args.planes),),
        grid_cache_path=Path(args.grid_cache),
        frequencies_mhz=frequencies_mhz,
    )
    problem = _subset_problem_cases(build_inverse_problem(config), [args.case_name])
    case = problem.cases[0]
    background_grid = problem.background_grids[0]
    grid = _apply_fit_params_to_grid(problem, background_grid, _clip_params(params))
    o_image = np.zeros((len(frequencies_mhz), problem.range_centers_km.size), dtype=float)
    x_image = np.zeros_like(o_image)
    doppler_num = np.zeros_like(o_image)
    doppler_den = np.zeros_like(o_image)
    for index, frequency_mhz in enumerate(frequencies_mhz):
        result = _render_frequency_plot_slice(problem, case, grid, float(frequency_mhz))
        o_image[index] = result.o_image
        x_image[index] = result.x_image
        doppler_num[index] = result.doppler_num
        doppler_den[index] = result.doppler_den
    np.savez_compressed(
        args.child_output,
        frequencies_mhz=np.asarray(frequencies_mhz, dtype=float),
        o_image=o_image,
        x_image=x_image,
        doppler_num=doppler_num,
        doppler_den=doppler_den,
        range_edges_km=np.asarray(problem.range_edges_km, dtype=float),
        range_centers_km=np.asarray(problem.range_centers_km, dtype=float),
    )


def _evaluate_params(
    *,
    script_path: Path,
    config: TopsideInverseConfig,
    case_name: str,
    params: IonosphereFitParams,
    block_size: int,
    workdir: Path,
) -> dict[str, np.ndarray]:
    frequencies = np.asarray(config.frequencies_mhz, dtype=float)
    nfreq = frequencies.size
    o_image = None
    x_image = None
    doppler_num = None
    doppler_den = None
    range_edges_km = None
    range_centers_km = None
    for start in range(0, nfreq, block_size):
        block = frequencies[start:start + block_size]
        output_path = workdir / f"block_{start:03d}.npz"
        command = [
            sys.executable,
            str(script_path),
            "--child-block",
            "--case-name",
            case_name,
            "--planes",
            str(config.planes[0]),
            "--grid-cache",
            str(config.grid_cache_path),
            "--child-frequencies-json",
            json.dumps(block.tolist()),
            "--child-params-json",
            json.dumps(vars(_clip_params(params))),
            "--child-output",
            str(output_path),
        ]
        subprocess.run(
            command,
            cwd=str(script_path.parent),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        payload = np.load(output_path)
        if o_image is None:
            nrange = int(payload["o_image"].shape[1])
            o_image = np.zeros((nfreq, nrange), dtype=float)
            x_image = np.zeros((nfreq, nrange), dtype=float)
            doppler_num = np.zeros((nfreq, nrange), dtype=float)
            doppler_den = np.zeros((nfreq, nrange), dtype=float)
            range_edges_km = np.asarray(payload["range_edges_km"], dtype=float)
            range_centers_km = np.asarray(payload["range_centers_km"], dtype=float)
        block_len = int(payload["o_image"].shape[0])
        o_image[start:start + block_len] = np.asarray(payload["o_image"], dtype=float)
        x_image[start:start + block_len] = np.asarray(payload["x_image"], dtype=float)
        doppler_num[start:start + block_len] = np.asarray(payload["doppler_num"], dtype=float)
        doppler_den[start:start + block_len] = np.asarray(payload["doppler_den"], dtype=float)
        output_path.unlink(missing_ok=True)
    return {
        "o_image": o_image,
        "x_image": x_image,
        "doppler_num": doppler_num,
        "doppler_den": doppler_den,
        "range_edges_km": range_edges_km,
        "range_centers_km": range_centers_km,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Memory-safe dense single-case fit driver.")
    parser.add_argument("--child-block", action="store_true")
    parser.add_argument("--child-frequencies-json", default=None)
    parser.add_argument("--child-params-json", default=None)
    parser.add_argument("--child-output", default=None)
    parser.add_argument("--case-name", default=DEFAULT_CASE)
    parser.add_argument("--planes", type=int, default=1)
    parser.add_argument("--grid-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--freq-start-mhz", type=float, default=1.0)
    parser.add_argument("--freq-stop-mhz", type=float, default=6.0)
    parser.add_argument("--freq-step-khz", type=float, default=100.0)
    parser.add_argument("--block-size", type=int, default=10)
    parser.add_argument("--fit-out", type=Path, default=Path("/tmp/plane01_sv121_oblique_dense_fit.json"))
    parser.add_argument("--screen-step-khz", type=float, default=1000.0)
    args = parser.parse_args()

    if args.child_block:
        _child_trace_block(args)
        return

    config = replace(
        TopsideInverseConfig(),
        planes=(int(args.planes),),
        grid_cache_path=args.grid_cache.expanduser(),
        frequencies_mhz=tuple(make_frequency_axis_mhz(args.freq_start_mhz, args.freq_stop_mhz, args.freq_step_khz).tolist()),
    )
    problem = _subset_problem_cases(build_inverse_problem(config), [args.case_name])
    case = problem.cases[0]
    los_bearing_deg, _, _ = _line_of_sight_angles(case.tx_points[1], case.rx_points[1])
    normal_bearings_deg = (_wrap_bearing_deg(los_bearing_deg - 90.0), _wrap_bearing_deg(los_bearing_deg + 90.0))

    script_path = Path(__file__).resolve()
    full_cache: dict[tuple[float, float, float, float, float], tuple[float, dict[str, np.ndarray]]] = {}
    screen_cache: dict[tuple[float, float, float, float, float], tuple[float, dict[str, np.ndarray]]] = {}

    with tempfile.TemporaryDirectory(prefix="dense_case_fit_") as temp_dir:
        workdir = Path(temp_dir)

        observed_raw = _evaluate_params(
            script_path=script_path,
            config=config,
            case_name=args.case_name,
            params=DEFAULT_TRUTH,
            block_size=max(int(args.block_size), 1),
            workdir=workdir,
        )
        observed_full = _finalize_raw(observed_raw, config)
        full_step_khz = max(float(args.freq_step_khz), 1e-6)
        screen_stride = max(1, int(round(float(args.screen_step_khz) / full_step_khz)))
        screen_indices = np.arange(0, len(config.frequencies_mhz), screen_stride, dtype=int)
        screen_config = replace(
            config,
            frequencies_mhz=tuple(np.asarray(config.frequencies_mhz, dtype=float)[screen_indices].tolist()),
        )
        observed_screen = _frequency_subset(observed_full, screen_indices)

        def evaluate(
            params: IonosphereFitParams,
            label: str,
            *,
            config_for_eval: TopsideInverseConfig,
            observed_for_eval: dict[str, np.ndarray],
            cache_for_eval: dict[tuple[float, float, float, float, float], tuple[float, dict[str, np.ndarray]]],
        ) -> tuple[float, dict[str, np.ndarray]]:
            clipped = _clip_params(params)
            key = _candidate_key(clipped)
            cached = cache_for_eval.get(key)
            if cached is not None:
                print(f"{label}: cost={cached[0]:.6f} params={vars(clipped)} cache=hit", flush=True)
                return cached
            raw = _evaluate_params(
                script_path=script_path,
                config=config_for_eval,
                case_name=args.case_name,
                params=clipped,
                block_size=max(int(args.block_size), 1),
                workdir=workdir,
            )
            finalized = _finalize_raw(raw, config_for_eval)
            cost = _single_case_cost(observed_for_eval, finalized)
            cache_for_eval[key] = (cost, finalized)
            print(f"{label}: cost={cost:.6f} params={vars(clipped)}", flush=True)
            return cost, finalized

        best_background = IonosphereFitParams(
            density_scale=1.0,
            hmf2_shift_km=0.0,
            wave_amplitude_fraction=0.0,
            wave_phase_rad=0.0,
            wave_bearing_deg=normal_bearings_deg[0],
        )
        best_background_cost, best_background_eval = evaluate(
            best_background,
            "background-seed",
            config_for_eval=config,
            observed_for_eval=observed_full,
            cache_for_eval=full_cache,
        )
        background_shift_candidates = (0.0, 15.0, 30.0)
        background_scale_candidates = (0.95, 1.05)
        for shift_km in background_shift_candidates[1:]:
            trial = replace(best_background, hmf2_shift_km=shift_km)
            cost, evaluated = evaluate(
                trial,
                "background-shift",
                config_for_eval=config,
                observed_for_eval=observed_full,
                cache_for_eval=full_cache,
            )
            if cost < best_background_cost:
                best_background_cost = cost
                best_background = _clip_params(trial)
                best_background_eval = evaluated
        for scale in background_scale_candidates:
            trial = replace(best_background, density_scale=scale)
            cost, evaluated = evaluate(
                trial,
                "background-scale",
                config_for_eval=config,
                observed_for_eval=observed_full,
                cache_for_eval=full_cache,
            )
            if cost < best_background_cost:
                best_background_cost = cost
                best_background = _clip_params(trial)
                best_background_eval = evaluated

        residual = np.asarray(observed_full["total_image"], dtype=float) - np.asarray(best_background_eval["total_image"], dtype=float)
        residual_power = float(np.sqrt(np.mean(residual * residual)))
        amplitude_seeds = (
            max(0.05, min(0.16, 1.0 * residual_power)),
            max(0.10, min(0.24, 1.8 * residual_power)),
        )
        amplitude_seeds = tuple(sorted({round(float(value), 4) for value in amplitude_seeds}))
        bearing_seeds = tuple(
            sorted(
                {
                    round(_wrap_bearing_deg(seed), 6)
                    for seed in (
                        normal_bearings_deg[0] - 45.0,
                        normal_bearings_deg[0],
                        normal_bearings_deg[0] + 45.0,
                        normal_bearings_deg[1] - 45.0,
                        normal_bearings_deg[1],
                        normal_bearings_deg[1] + 45.0,
                    )
                }
            )
        )
        phase_seeds = (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)

        screen_ranked: list[tuple[float, IonosphereFitParams]] = []
        seen_wave: set[tuple[float, float, float, float, float]] = set()
        for amplitude in amplitude_seeds:
            for bearing_deg in bearing_seeds:
                for phase_rad in phase_seeds:
                    trial = replace(
                        best_background,
                        wave_amplitude_fraction=amplitude,
                        wave_phase_rad=phase_rad,
                        wave_bearing_deg=bearing_deg,
                    )
                    key = _candidate_key(trial)
                    if key in seen_wave:
                        continue
                    seen_wave.add(key)
                    cost, _ = evaluate(
                        trial,
                        "wave-screen",
                        config_for_eval=screen_config,
                        observed_for_eval=observed_screen,
                        cache_for_eval=screen_cache,
                    )
                    screen_ranked.append((cost, _clip_params(trial)))
        screen_ranked.sort(key=lambda item: item[0])
        full_wave_ranked: list[tuple[float, IonosphereFitParams]] = []
        for _, candidate in screen_ranked[:2]:
            cost, _ = evaluate(
                candidate,
                "wave-dense",
                config_for_eval=config,
                observed_for_eval=observed_full,
                cache_for_eval=full_cache,
            )
            full_wave_ranked.append((cost, candidate))
        best_wave_cost, best_wave = min(full_wave_ranked, key=lambda item: item[0])
        for trial in (
            replace(best_wave, density_scale=best_wave.density_scale - 0.05),
            replace(best_wave, density_scale=best_wave.density_scale + 0.05),
            replace(best_wave, hmf2_shift_km=best_wave.hmf2_shift_km - 10.0),
            replace(best_wave, hmf2_shift_km=best_wave.hmf2_shift_km + 10.0),
            replace(best_wave, wave_amplitude_fraction=max(0.04, best_wave.wave_amplitude_fraction - 0.04)),
            replace(best_wave, wave_amplitude_fraction=min(0.25, best_wave.wave_amplitude_fraction + 0.04)),
            replace(best_wave, wave_phase_rad=best_wave.wave_phase_rad - 2.0 * math.pi / 3.0),
            replace(best_wave, wave_phase_rad=best_wave.wave_phase_rad + 2.0 * math.pi / 3.0),
            replace(best_wave, wave_bearing_deg=best_wave.wave_bearing_deg - 15.0),
            replace(best_wave, wave_bearing_deg=best_wave.wave_bearing_deg + 15.0),
        ):
            cost, _ = evaluate(
                trial,
                "wave-refine",
                config_for_eval=config,
                observed_for_eval=observed_full,
                cache_for_eval=full_cache,
            )
            if cost < best_wave_cost:
                best_wave_cost = cost
                best_wave = _clip_params(trial)

        payload = {
            "case_name": args.case_name,
            "frequencies_mhz": list(map(float, config.frequencies_mhz)),
            "truth_params": vars(DEFAULT_TRUTH),
            "start_guess_method": {
                "background": "dense 100 kHz background-only search over hmf2_shift_km and density_scale using the observed dense ionogram only",
                "wave": "screen nonzero-wave seeds against a downsampled observed ionogram, then reevaluate the best seeds on the full dense ionogram",
            },
            "initial_params": vars(IonosphereFitParams(density_scale=1.0, hmf2_shift_km=0.0, wave_amplitude_fraction=0.0, wave_phase_rad=0.0, wave_bearing_deg=normal_bearings_deg[0])),
            "screen_step_khz": float(args.screen_step_khz),
            "residual_rms_after_background": residual_power,
            "best_background_only": _candidate_record(best_background_cost, best_background),
            "best_wave_constrained": _candidate_record(best_wave_cost, best_wave),
            "wave_screen_top2": [_candidate_record(cost, params) for cost, params in screen_ranked[:2]],
            "full_dense_evaluations": len(full_cache),
            "screen_evaluations": len(screen_cache),
            "normal_bearings_deg": [float(value) for value in normal_bearings_deg],
        }
        args.fit_out.expanduser().write_text(json.dumps(payload, indent=2))
        print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
