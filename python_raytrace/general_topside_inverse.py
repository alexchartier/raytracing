"""Fast vertical topside inversion with an explicit, flexible F2 peak.

The 1-D group-range operator is a search surrogate. Finalists intended for a
three-dimensional retrieval must be checked with the full O/X ray tracer.
No truth density enters the fit; this module accepts only ionogram returns.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
from functools import cached_property

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, least_squares


PLASMA_MHZ_PER_SQRT_CM3 = 0.00898
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)


def vertical_group_range_km(profile: object, frequency_mhz: float,
                            spacecraft_altitude_km: float) -> float:
    """Two-way vertical scalar-plasma group range to the upper F2 turning point."""
    f = float(frequency_mhz)
    z_sc = float(spacecraft_altitude_km)
    if not profile.hmf2_km < z_sc or not 0 < f < profile.fof2_mhz:
        return float("nan")
    f_sc = PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(profile.density_cm3(z_sc))
    if f <= f_sc:
        return float("nan")

    def q(z: np.ndarray | float) -> np.ndarray:
        return 1.0 - (PLASMA_MHZ_PER_SQRT_CM3 ** 2 / f**2) * profile.density_cm3(z)

    turning = brentq(lambda z: float(q(z)), profile.hmf2_km, z_sc,
                     xtol=1e-8)
    # z = turning + u^2 removes the square-root singularity at reflection.
    umax = np.sqrt(z_sc - turning)
    u = 0.5 * umax * (_NODES + 1.0)
    integrand = 2.0 * u / np.sqrt(np.maximum(q(turning + u * u), 1e-14))
    return float(umax * np.dot(_WEIGHTS, integrand))


@dataclass(frozen=True)
class ChapmanF2:
    fof2_mhz: float
    hmf2_km: float
    bottomside_scale_km: float
    topside_scale_km: float
    tail_curvature: float = 0.0

    def density_cm3(self, altitudes_km: np.ndarray | float) -> np.ndarray:
        """Generalized Chapman layer; curvature changes only the upper tail."""
        z = np.asarray(altitudes_km, dtype=float)
        dz = z - self.hmf2_km
        below = dz < 0.0
        t = np.where(below, dz / self.bottomside_scale_km,
                     dz / self.topside_scale_km)
        top = np.maximum(t, 0.0)
        # This transform preserves the peak location, peak curvature, and a
        # monotone topside for tail_curvature > -1.
        denominator = 1.0 + self.tail_curvature * top / (top + 3.0)
        y = np.where(below, t, t / denominator)
        exponent = 0.5 * (1.0 - y - np.exp(np.clip(-y, -80.0, 80.0)))
        return (self.fof2_mhz / PLASMA_MHZ_PER_SQRT_CM3) ** 2 * np.exp(exponent)

    def group_range_km(self, frequency_mhz: float, spacecraft_altitude_km: float) -> float:
        return vertical_group_range_km(self, frequency_mhz, spacecraft_altitude_km)


@dataclass(frozen=True)
class IRITopsideBasis:
    """Peak-aligned log-density PCA from geographically varied IRI profiles."""

    offset_km: np.ndarray
    mean_log_ratio: np.ndarray
    components: np.ndarray  # Rows have one sample-standard-deviation amplitude.
    source: str
    profile_count: int

    @classmethod
    def read(cls, path: Path) -> "IRITopsideBasis":
        with np.load(path, allow_pickle=False) as data:
            return cls(data["offset_km"], data["mean_log_ratio"],
                       data["components"], str(data["source"]),
                       int(data["profile_count"]))

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, offset_km=self.offset_km,
                            mean_log_ratio=self.mean_log_ratio,
                            components=self.components, source=np.array(self.source),
                            profile_count=np.array(self.profile_count))

    def profile(self, fof2_mhz: float, hmf2_km: float,
                scores: np.ndarray, bottomside_scale_km: float = 55.0) -> "IRIBasisF2":
        return IRIBasisF2(float(fof2_mhz), float(hmf2_km),
                          float(bottomside_scale_km), self,
                          np.asarray(scores, dtype=float))


@dataclass(frozen=True)
class IRIBasisF2:
    fof2_mhz: float
    hmf2_km: float
    bottomside_scale_km: float
    basis: IRITopsideBasis
    scores: np.ndarray

    @cached_property
    def _interpolator(self) -> PchipInterpolator:
        shape = self.basis.mean_log_ratio + self.scores @ self.basis.components
        # Physical topside profiles fall from the F2 maximum. The optimizer
        # penalizes violations, and this projection keeps trial ray integrals
        # single-turning-point even at a bound or during finite differencing.
        shape = np.minimum.accumulate(shape)
        return PchipInterpolator(self.basis.offset_km, shape, extrapolate=True)

    def density_cm3(self, altitudes_km: np.ndarray | float) -> np.ndarray:
        z = np.asarray(altitudes_km, dtype=float)
        offset = z - self.hmf2_km
        top = np.exp(np.minimum(self._interpolator(np.maximum(offset, 0.0)), 0.0))
        t = offset / self.bottomside_scale_km
        bottom = np.exp(0.5 * (1.0 - t - np.exp(np.clip(-t, -80.0, 80.0))))
        return (self.fof2_mhz / PLASMA_MHZ_PER_SQRT_CM3) ** 2 * np.where(offset >= 0, top, bottom)

    def group_range_km(self, frequency_mhz: float, spacecraft_altitude_km: float) -> float:
        return vertical_group_range_km(self, frequency_mhz, spacecraft_altitude_km)


SPLINE_OFFSETS_KM = np.array([0.0, 20.0, 50.0, 100.0, 180.0,
                              280.0, 400.0, 550.0, 700.0])


@dataclass(frozen=True)
class MonotoneSplineF2:
    """Model-agnostic topside log density with positive interval slopes."""

    fof2_mhz: float
    hmf2_km: float
    bottomside_scale_km: float
    log_slope_per_km: np.ndarray

    @cached_property
    def _interpolator(self) -> PchipInterpolator:
        slope = np.asarray(self.log_slope_per_km, dtype=float)
        if slope.shape != (len(SPLINE_OFFSETS_KM) - 1,) or np.any(slope <= 0):
            raise ValueError("Expected positive topside interval slopes")
        log_ratio = np.r_[0.0, -np.cumsum(slope * np.diff(SPLINE_OFFSETS_KM))]
        return PchipInterpolator(SPLINE_OFFSETS_KM, log_ratio, extrapolate=True)

    def density_cm3(self, altitudes_km: np.ndarray | float) -> np.ndarray:
        z = np.asarray(altitudes_km, dtype=float)
        offset = z - self.hmf2_km
        top = np.exp(np.minimum(self._interpolator(np.maximum(offset, 0.0)), 0.0))
        t = offset / self.bottomside_scale_km
        bottom = np.exp(0.5 * (1.0 - t - np.exp(np.clip(-t, -80.0, 80.0))))
        return (self.fof2_mhz / PLASMA_MHZ_PER_SQRT_CM3) ** 2 * np.where(offset >= 0, top, bottom)

    def group_range_km(self, frequency_mhz: float, spacecraft_altitude_km: float) -> float:
        return vertical_group_range_km(self, frequency_mhz, spacecraft_altitude_km)


def build_global_iri2020_basis() -> IRITopsideBasis:
    """Generate an IRI-2020 prior without accessing any validation truth grid."""
    from python_raytrace.iri2020_model import Iri2020Bridge, r12_from_f107

    model = Iri2020Bridge()
    altitude = np.arange(100.0, 1200.1, 5.0)
    offset = np.arange(0.0, 600.1, 5.0)
    samples = []
    # Latitude, longitude, season, and local time all vary. The exact D-case
    # position and epoch are not among these training points.
    epochs = (datetime(2011, 3, 20, 0), datetime(2011, 6, 21, 6),
              datetime(2011, 9, 22, 12), datetime(2011, 12, 21, 18))
    for latitude in (-70.0, -35.0, 0.0, 35.0, 70.0):
        for longitude in (-144.0, -72.0, 0.0, 72.0, 144.0):
            for epoch in epochs:
                result = model.profile(epoch, latitude, longitude,
                                       r12=r12_from_f107(120.0),
                                       altitudes_km=altitude, d_region_model="none")
                density = result.electron_density_m3 / 1e6
                if not np.all(np.isfinite(density)) or np.any(density <= 0):
                    continue
                maximum = int(np.argmax(density))
                peak_height = altitude[maximum]
                if peak_height < 180.0 or peak_height > 550.0:
                    continue
                log_density = np.log(density)
                curve = np.interp(peak_height + offset, altitude, log_density)
                curve -= curve[0]
                if np.any(np.diff(curve) > 0.025):
                    continue
                samples.append(curve)
    if len(samples) < 50:
        raise RuntimeError("IRI basis sampling produced too few valid profiles")
    matrix = np.asarray(samples)
    mean = np.mean(matrix, axis=0)
    _, singular, vt = np.linalg.svd(matrix - mean, full_matrices=False)
    components = vt[:4] * (singular[:4, None] / np.sqrt(len(samples) - 1))
    return IRITopsideBasis(offset, mean, components,
                           "IRI-2020 bridge; 2011 global 5x5x4 sample; F10.7=120", len(samples))


@dataclass(frozen=True)
class FitResult:
    layer: ChapmanF2 | IRIBasisF2 | MonotoneSplineF2
    x_frequency_offset_mhz: float
    x_range_offset_km: float
    objective: float
    range_mae_km: dict[int, float]
    fitted_frequencies: dict[int, int]


def fit_spline_ionogram(
    records: np.ndarray,
    frequencies_mhz: np.ndarray,
    spacecraft_altitude_km: float,
    *,
    bottomside_scale_km: float = 55.0,
    regularize_tail: bool = False,
    fof2_anchor_mhz: float | None = None,
) -> FitResult:
    """Fit O-mode group ranges with a free, monotone topside profile.

    The shape has eight independent log-density slopes and no empirical-model
    profile in its basis. With ``regularize_tail``, positive upper gradients
    and smooth log slopes keep the unobserved tail from becoming a plateau.
    Full O/X rays must score the fitted profile.
    """
    freqs = np.asarray(frequencies_mhz, dtype=float)
    observed_f, observed_r = _ridges(records, freqs, 1)
    if len(observed_f) < 10:
        raise ValueError("At least ten O-mode frequency bins are required")
    df = float(np.median(np.diff(freqs)))
    nose = float(observed_f[-1])
    anchor = nose + 0.5 * df if fof2_anchor_mhz is None else float(fof2_anchor_mhz)
    if not np.isfinite(anchor) or anchor < nose - df:
        raise ValueError("foF2 anchor must be finite and near or above the O nose")
    keep = observed_f <= nose - df
    fit_f, fit_r = observed_f[keep], observed_r[keep]
    if len(fit_f) < 8:
        raise ValueError("Too few O-mode bins below the nose")
    offsets = SPLINE_OFFSETS_KM

    def initial_slopes(scale_km: float) -> np.ndarray:
        seed = ChapmanF2(nose + 0.05, 270.0, bottomside_scale_km, scale_km)
        values = seed.density_cm3(270.0 + offsets)
        return -np.diff(np.log(values)) / np.diff(offsets)

    slope_lower = np.full(len(offsets) - 1, 0.00002)
    if regularize_tail:
        # The sweep begins at 2 MHz: it cannot resolve plasma much above the
        # corresponding turning height. Enforce only a very broad, monotone
        # decay in the extrapolated tail (e-folding length <= 500 km).
        slope_lower[-3:] = 0.002
    lower_fo = max(1.0, nose - 0.05, anchor - 0.25)
    upper_fo = min(15.0, max(nose + 0.50, anchor + 0.25))
    lower = np.r_[lower_fo, 170.0, slope_lower]
    upper = np.r_[upper_fo,
                  min(500.0, spacecraft_altitude_km - 60.0),
                  np.full(len(offsets) - 1, 0.05)]

    def make_profile(x: np.ndarray) -> MonotoneSplineF2:
        return MonotoneSplineF2(float(x[0]), float(x[1]),
                                bottomside_scale_km, np.asarray(x[2:]))

    def residual(x: np.ndarray) -> np.ndarray:
        layer = make_profile(x)
        predicted = np.array([layer.group_range_km(f, spacecraft_altitude_km)
                              for f in fit_f])
        invalid = ~np.isfinite(predicted)
        predicted[invalid] = fit_r[invalid] + 500.0
        slopes = x[2:]
        # Normalize adjacent slope differences by their characteristic
        # topside gradient; this remains weak relative to a 30-km range miss.
        smoothness = np.diff(slopes) / 0.018
        return np.r_[(predicted - fit_r) / 30.0,
                     (x[0] - anchor) / (0.20 if fof2_anchor_mhz is None else 0.15),
                     0.35 * smoothness]

    trials = []
    for scale in (55.0, 85.0, 125.0):
        start = np.r_[anchor, 270.0, initial_slopes(scale)]
        start = np.clip(start, lower + 1e-6, upper - 1e-6)
        fitted = least_squares(residual, start, bounds=(lower, upper),
                               loss="soft_l1", f_scale=1.0, max_nfev=180,
                               diff_step=1e-4)
        trials.append(fitted)
    best = min(trials, key=lambda result: float(np.mean(residual(result.x)**2)))
    profile = make_profile(best.x)
    predicted = np.array([profile.group_range_km(f, spacecraft_altitude_km)
                          for f in fit_f])
    return FitResult(profile, 0.0, 0.0,
                     float(np.mean(residual(best.x)**2)),
                     {1: float(np.mean(np.abs(predicted - fit_r)))},
                     {1: len(fit_f)})


def fit_iri_basis_ionogram(
    records: np.ndarray,
    frequencies_mhz: np.ndarray,
    spacecraft_altitude_km: float,
    basis: IRITopsideBasis,
    *,
    fit_x_mode: bool = True,
    initial_hmf2_km: float = 280.0,
) -> FitResult:
    """Fit explicit foF2/hmF2 and four IRI-derived topside shape scores."""
    freqs = np.asarray(frequencies_mhz, dtype=float)
    ridges = {m: _ridges(records, freqs, m) for m in (1, -1)}
    of, _ = ridges[1]
    if len(of) < 8:
        raise ValueError("At least eight O-mode frequency bins are required")
    nose = float(of[-1])
    df = float(np.median(np.diff(freqs)))
    nshape = len(basis.components)
    x0 = np.r_[nose + df, initial_hmf2_km, np.zeros(nshape), 0.15, 0.0]
    lower = np.r_[nose + 0.001, 170.0, np.full(nshape, -2.5), -0.1, -100.0]
    upper = np.r_[min(15.0, nose + 0.45), min(500.0, spacecraft_altitude_km - 60.0),
                   np.full(nshape, 2.5), 0.6, 100.0]
    x0 = np.clip(x0, lower + 1e-4, upper - 1e-4)
    data = {}
    for mode in (1, -1) if fit_x_mode else (1,):
        f, r = ridges[mode]
        keep = f <= min(float(f[-1]) - df, nose - df)
        data[mode] = (f[keep], r[keep])

    def make_profile(x: np.ndarray) -> IRIBasisF2:
        return basis.profile(x[0], x[1], x[2:2+nshape])

    def residual(x: np.ndarray) -> np.ndarray:
        layer = make_profile(x)
        out = []
        for mode, (f, observed) in data.items():
            effective = f - x[-2] if mode == -1 else f
            predicted = np.array([layer.group_range_km(ff, spacecraft_altitude_km)
                                  for ff in effective])
            if mode == -1:
                predicted += x[-1]
            invalid = ~np.isfinite(predicted)
            predicted[invalid] = observed[invalid] + 500.0
            out.extend((predicted - observed) / 30.0)
        # One sample past a 100-kHz nose is a plausible critical frequency.
        # Keep it soft because a homing gap can hide real high-frequency rays.
        out.append((x[0] - (nose + 0.5 * df)) / 0.20)
        out.extend(x[2:2+nshape] / 1.5)
        shape = layer.basis.mean_log_ratio + x[2:2+nshape] @ layer.basis.components
        out.extend(np.maximum(np.diff(shape), 0.0) / 0.02)
        out.extend([x[-2] / 0.25, x[-1] / 50.0])
        return np.asarray(out)

    fitted = least_squares(residual, x0, bounds=(lower, upper), loss="soft_l1",
                           f_scale=1.0, max_nfev=200, diff_step=1e-4)
    profile = make_profile(fitted.x)
    mae = {}
    counts = {}
    for mode, (f, observed) in data.items():
        effective = f - fitted.x[-2] if mode == -1 else f
        predicted = np.array([profile.group_range_km(ff, spacecraft_altitude_km)
                              for ff in effective])
        if mode == -1:
            predicted += fitted.x[-1]
        mae[mode] = float(np.mean(np.abs(predicted - observed)))
        counts[mode] = len(f)
    return FitResult(profile, float(fitted.x[-2]), float(fitted.x[-1]),
                     float(np.mean(residual(fitted.x) ** 2)), mae, counts)


def _ridges(records: np.ndarray, frequencies: np.ndarray, mode: int) -> tuple[np.ndarray, np.ndarray]:
    rows = np.asarray(records, dtype=float)
    rows = rows[rows[:, 1] == mode]
    indices = np.unique(rows[:, 0]).astype(int)
    ranges = np.array([np.median(rows[rows[:, 0] == index, 2]) for index in indices])
    return frequencies[indices], ranges


def fit_vertical_ionogram(
    records: np.ndarray,
    frequencies_mhz: np.ndarray,
    spacecraft_altitude_km: float,
    *,
    initial_fof2_mhz: float | None = None,
    bottomside_scale_km: float = 55.0,
    fit_x_mode: bool = True,
) -> FitResult:
    """Fit O/X ridges with a robust loss and a weak prior on tail curvature.

    The X operator is a nuisance frequency and range correction to the O
    surrogate, since a one-dimensional scalar plasma model omits magnetoionic
    refraction. These corrections must not be interpreted as plasma parameters.
    Frequencies above the accepted-return nose are treated as censored, not as
    proof that no reflection exists; the 100-kHz nose has finite resolution.
    """
    freqs = np.asarray(frequencies_mhz, dtype=float)
    if freqs.ndim != 1 or len(freqs) < 3 or not np.all(np.diff(freqs) > 0):
        raise ValueError("frequencies_mhz must be a sorted frequency axis")
    ridges = {mode: _ridges(records, freqs, mode) for mode in (1, -1)}
    of, _ = ridges[1]
    if len(of) < 8:
        raise ValueError("At least eight O-mode frequency bins are required")
    nose = float(of[-1])
    df = float(np.median(np.diff(freqs)))
    fof_init = initial_fof2_mhz if initial_fof2_mhz is not None else nose + 0.5 * df
    x0 = np.array([fof_init, 270.0, 75.0, 0.0, 0.15, 0.0], dtype=float)
    lower = np.array([max(1.0, nose - 0.05), 170.0, 20.0, -0.75, -0.1, -100.0])
    upper = np.array([min(15.0, nose + 0.45), min(500.0, spacecraft_altitude_km - 60.0),
                      220.0, 2.0, 0.6, 100.0])
    x0 = np.clip(x0, lower + 1e-4, upper - 1e-4)

    def unpack(x: np.ndarray) -> ChapmanF2:
        return ChapmanF2(float(x[0]), float(x[1]), bottomside_scale_km,
                         float(x[2]), float(x[3]))

    # Use all accepted frequency bins except the last two, where singular
    # group delay and occasional failed homing make the 1-D proxy least stable.
    data = {}
    for mode in (1, -1) if fit_x_mode else (1,):
        f, r = ridges[mode]
        keep = f <= min(float(f[-1]) - 2 * df, nose - df)
        data[mode] = (f[keep], r[keep])

    def residual(x: np.ndarray) -> np.ndarray:
        layer = unpack(x)
        out = []
        for mode, (f, observed) in data.items():
            effective = f - x[4] if mode == -1 else f
            predicted = np.array([layer.group_range_km(ff, spacecraft_altitude_km)
                                  for ff in effective])
            if mode == -1:
                predicted += x[5]
            invalid = ~np.isfinite(predicted)
            predicted[invalid] = observed[invalid] + 500.0
            # The full-ray median may have several branches. Huber-like loss
            # at 30 km prevents one mismatched branch from setting hmF2.
            out.extend((predicted - observed) / 30.0)
        # Nose supplies a bounded foF2 constraint; X nuisance quantities and
        # curvature remain regularized unless returns genuinely determine them.
        out.extend([(x[0] - (nose + 0.5 * df)) / 0.20,
                    x[3] / 0.75, x[4] / 0.25, x[5] / 50.0])
        return np.asarray(out)

    fit = least_squares(residual, x0, bounds=(lower, upper), loss="soft_l1",
                        f_scale=1.0, max_nfev=170, diff_step=1e-4)
    layer = unpack(fit.x)
    mae = {}
    counts = {}
    for mode, (f, observed) in data.items():
        effective = f - fit.x[4] if mode == -1 else f
        predicted = np.array([layer.group_range_km(ff, spacecraft_altitude_km)
                              for ff in effective])
        if mode == -1:
            predicted += fit.x[5]
        mae[mode] = float(np.mean(np.abs(predicted - observed)))
        counts[mode] = len(f)
    return FitResult(layer, float(fit.x[4]), float(fit.x[5]),
                     float(np.mean(residual(fit.x) ** 2)), mae, counts)
