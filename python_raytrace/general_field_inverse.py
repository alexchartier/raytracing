"""Model-independent density-field corrections for full-ray topside inversion.

The prior may come from any ionosphere model. A small, smooth correction to
log electron density provides peak, height, and topside-shape freedom along a
spacecraft pass. A measured local density is imposed separately in the upper
tail. Candidate selection must use forward-traced O/X ionograms.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace

import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator

from .absorption import effective_collision_frequency
from .grid import IonosphereGrid

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class LocalDensity:
    """In-situ observations along one pass at a common spacecraft altitude."""

    lat_lon_deg: np.ndarray  # (sample, 2)
    altitude_km: float
    electron_density_cm3: np.ndarray
    relative_uncertainty: float = 0.005

    def __post_init__(self) -> None:
        locations = np.asarray(self.lat_lon_deg, dtype=float)
        density = np.asarray(self.electron_density_cm3, dtype=float)
        if (locations.ndim != 2 or locations.shape[1] != 2
                or density.shape != (len(locations),) or not len(locations)
                or not np.all(np.isfinite(locations))
                or not np.all(np.isfinite(density)) or np.any(density <= 0)
                or not np.isfinite(self.altitude_km)
                or not np.isfinite(self.relative_uncertainty)
                or not 0 < self.relative_uncertainty < 1):
            raise ValueError("Expected finite locations and positive local densities")


@dataclass(frozen=True)
class TrackFrame:
    latitude_origin_deg: float
    longitude_origin_deg: float
    east_axis: float
    north_axis: float
    sample_distance_km: np.ndarray

    @classmethod
    def from_locations(cls, locations: np.ndarray) -> "TrackFrame":
        points = np.asarray(locations, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or not len(points):
            raise ValueError("Expected (sample, latitude/longitude) coordinates")
        lat0 = float(np.mean(points[:, 0]))
        lon0 = float(np.mean(np.rad2deg(np.unwrap(np.deg2rad(points[:, 1])))))
        east = EARTH_RADIUS_KM * np.cos(np.deg2rad(lat0)) * np.deg2rad(
            (points[:, 1] - lon0 + 180.0) % 360.0 - 180.0)
        north = EARTH_RADIUS_KM * np.deg2rad(points[:, 0] - lat0)
        xy = np.column_stack((east, north))
        if len(points) == 1 or np.max(np.linalg.norm(xy - xy[0], axis=1)) < 1e-6:
            axis = np.array([0.0, 1.0])
        else:
            _, _, vt = np.linalg.svd(xy - np.mean(xy, axis=0), full_matrices=False)
            axis = vt[0]
            if (xy[-1] - xy[0]) @ axis < 0:
                axis = -axis
        distance = xy @ axis
        return cls(lat0, lon0, float(axis[0]), float(axis[1]), distance)

    def project(self, latitude_deg: np.ndarray, longitude_deg: np.ndarray) -> np.ndarray:
        east = EARTH_RADIUS_KM * np.cos(np.deg2rad(self.latitude_origin_deg)) * np.deg2rad(
            (np.asarray(longitude_deg) - self.longitude_origin_deg + 180.0) % 360.0 - 180.0)
        north = EARTH_RADIUS_KM * np.deg2rad(
            np.asarray(latitude_deg) - self.latitude_origin_deg)
        return self.east_axis * east + self.north_axis * north


@dataclass(frozen=True)
class TrackLogDensityBasis:
    """Smooth along-track × peak-relative-altitude log-density basis."""

    frame: TrackFrame
    center_distance_km: np.ndarray
    spatial_width_km: float
    altitude_offsets_km: np.ndarray
    altitude_widths_km: np.ndarray

    @classmethod
    def from_locations(cls, locations: np.ndarray, *, spatial_centers: int = 4,
                       altitude_offsets_km: tuple[float, ...] =
                       (-50.0, 0.0, 80.0, 210.0, 420.0),
                       altitude_widths_km: tuple[float, ...] =
                       (55.0, 55.0, 75.0, 120.0, 190.0)) -> "TrackLogDensityBasis":
        frame = TrackFrame.from_locations(locations)
        s = frame.sample_distance_km
        if spatial_centers < 1 or len(altitude_offsets_km) != len(altitude_widths_km):
            raise ValueError("Invalid basis dimensions")
        if spatial_centers > 1 and np.ptp(s) < 1e-3:
            raise ValueError("Multiple spatial centers need a nonzero track span")
        centers = np.linspace(np.min(s), np.max(s), spatial_centers)
        width = (max(80.0, 1.2 * (centers[1] - centers[0]))
                 if spatial_centers > 1 else 500.0)
        return cls(frame, centers, float(width),
                   np.asarray(altitude_offsets_km, dtype=float),
                   np.asarray(altitude_widths_km, dtype=float))

    @property
    def coefficient_shape(self) -> tuple[int, int]:
        return len(self.center_distance_km), len(self.altitude_offsets_km)

    def correction(self, grid: IonosphereGrid, coefficients: np.ndarray) -> np.ndarray:
        c = np.asarray(coefficients, dtype=float)
        if c.shape != self.coefficient_shape or not np.all(np.isfinite(c)):
            raise ValueError(f"Expected finite coefficients shaped {self.coefficient_shape}")
        lat, lon = np.meshgrid(grid.latitudes_deg, grid.longitudes_deg, indexing="ij")
        distance = self.frame.project(lat, lon)
        log_spatial = -0.5 * ((distance[..., None] - self.center_distance_km)
                              / self.spatial_width_km) ** 2
        spatial = np.exp(log_spatial - np.max(log_spatial, axis=-1, keepdims=True))
        spatial /= np.sum(spatial, axis=-1, keepdims=True)
        peak = np.asarray(grid.altitudes_km)[np.argmax(grid.iono_en_grid, axis=2)]
        relative_altitude = (np.asarray(grid.altitudes_km)[None, None, :, None]
                             - peak[:, :, None, None])
        vertical = np.exp(-0.5 * ((relative_altitude - self.altitude_offsets_km)
                                  / self.altitude_widths_km) ** 2)
        return np.einsum("ija,ab,ijkb->ijk", spatial, c, vertical, optimize=True)


def sampled_local_density(grid: IonosphereGrid, measurement: LocalDensity) -> np.ndarray:
    points = np.column_stack((measurement.lat_lon_deg,
                              np.full(len(measurement.lat_lon_deg),
                                      measurement.altitude_km)))
    interpolation = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        grid.iono_en_grid, bounds_error=True)
    return np.asarray(interpolation(points), dtype=float)


def _anchor_upper_tail(grid: IonosphereGrid, density: np.ndarray,
                       measurement: LocalDensity, frame: TrackFrame,
                       *, taper_depth_km: float = 200.0,
                       relative_tolerance: float = 0.005) -> np.ndarray:
    """Interpolate local log-density residuals along track and taper below SC."""
    if taper_depth_km <= 0 or measurement.altitude_km <= grid.altitudes_km[0]:
        raise ValueError("Invalid local-density taper")
    if measurement.altitude_km > grid.altitudes_km[-1]:
        raise ValueError("Spacecraft lies above the candidate grid")
    sample_distance = frame.sample_distance_km
    order = np.argsort(sample_distance)
    sorted_distance = sample_distance[order]
    if len(sorted_distance) > 1 and np.min(np.diff(sorted_distance)) < 1e-4:
        raise ValueError("Local observations need distinct along-track positions")
    lat, lon = np.meshgrid(grid.latitudes_deg, grid.longitudes_deg, indexing="ij")
    grid_distance = frame.project(lat, lon)
    start = measurement.altitude_km - taper_depth_km
    fraction = np.clip((np.asarray(grid.altitudes_km) - start) / taper_depth_km,
                       0.0, 1.0)
    vertical_weight = fraction * fraction * (3.0 - 2.0 * fraction)
    corrected = density.copy()
    for _ in range(4):
        trial = replace(grid, iono_en_grid=corrected)
        modeled = sampled_local_density(trial, measurement)
        residual = np.log(measurement.electron_density_cm3 / modeled)
        if np.max(abs(np.expm1(residual))) <= relative_tolerance:
            return corrected
        if len(residual) == 1:
            along_track = np.full_like(grid_distance, residual[0], dtype=float)
        else:
            interpolator = PchipInterpolator(sorted_distance, residual[order])
            along_track = interpolator(np.clip(grid_distance, sorted_distance[0],
                                                sorted_distance[-1]))
        corrected *= np.exp(along_track[..., None] * vertical_weight[None, None, :])
    # The smooth interpolation can alias a wave onto a coarse forward grid.
    # Solve the remaining constraint at grid nodes. Prefer a whole latitude or
    # longitude row, preserving smoothness across the perpendicular direction.
    def bracketing(axis: np.ndarray, value: float) -> tuple[int, int, float]:
        upper = int(np.searchsorted(axis, value, side="right"))
        upper = min(max(upper, 1), len(axis) - 1)
        lower = upper - 1
        fraction = (value - axis[lower]) / (axis[upper] - axis[lower])
        if not 0.0 <= fraction <= 1.0:
            raise ValueError("Local observation falls outside the density grid")
        return lower, upper, float(fraction)

    lat_parts = [bracketing(np.asarray(grid.latitudes_deg), float(lat))
                 for lat in measurement.lat_lon_deg[:, 0]]
    lon_parts = [bracketing(np.asarray(grid.longitudes_deg), float(lon))
                 for lon in measurement.lat_lon_deg[:, 1]]
    span = np.ptp(measurement.lat_lon_deg, axis=0)
    prefer_latitude = (abs(frame.north_axis) >= abs(frame.east_axis)
                       and span[1] < np.median(np.diff(grid.longitudes_deg)))
    prefer_longitude = (abs(frame.east_axis) > abs(frame.north_axis)
                        and span[0] < np.median(np.diff(grid.latitudes_deg)))

    for axis_name in (["latitude", "nodes"] if prefer_latitude else
                      ["longitude", "nodes"] if prefer_longitude else ["nodes"]):
        trial_density = corrected.copy()
        at_spacecraft = np.empty(corrected.shape[:2])
        for i in range(corrected.shape[0]):
            for j in range(corrected.shape[1]):
                at_spacecraft[i, j] = np.interp(measurement.altitude_km,
                                                grid.altitudes_km, corrected[i, j])
        sample_nodes: list[list[tuple[int, float]]] = []
        for (i0, i1, ti), (j0, j1, tj) in zip(lat_parts, lon_parts):
            sample_nodes.append([
                (i0 * corrected.shape[1] + j0, (1 - ti) * (1 - tj)),
                (i0 * corrected.shape[1] + j1, (1 - ti) * tj),
                (i1 * corrected.shape[1] + j0, ti * (1 - tj)),
                (i1 * corrected.shape[1] + j1, ti * tj),
            ])
        if axis_name == "latitude":
            columns = list(range(corrected.shape[0]))
            node_column = lambda index: index // corrected.shape[1]
        elif axis_name == "longitude":
            columns = list(range(corrected.shape[1]))
            node_column = lambda index: index % corrected.shape[1]
        else:
            columns = sorted({index for row in sample_nodes for index, _ in row})
            node_column = lambda index: index
        column_index = {value: index for index, value in enumerate(columns)}
        matrix = np.zeros((len(sample_nodes), len(columns)))
        flat_density = at_spacecraft.ravel()
        for row_index, nodes in enumerate(sample_nodes):
            for node, weight in nodes:
                matrix[row_index, column_index[node_column(node)]] += (
                    weight * flat_density[node])
        delta, *_ = np.linalg.lstsq(
            matrix, measurement.electron_density_cm3 - matrix @ np.ones(len(columns)),
            rcond=1e-10)
        factor = 1.0 + delta
        if np.any(factor <= 0):
            continue
        if axis_name == "latitude":
            log_factor = np.log(factor)[:, None]
        elif axis_name == "longitude":
            log_factor = np.log(factor)[None, :]
        else:
            log_factor = np.zeros(corrected.shape[:2])
            for node, value in zip(columns, np.log(factor)):
                log_factor[np.unravel_index(node, log_factor.shape)] = value
        trial_density *= np.exp(log_factor[..., None] * vertical_weight[None, None, :])
        remaining = (sampled_local_density(replace(grid, iono_en_grid=trial_density),
                                           measurement)
                     / measurement.electron_density_cm3 - 1)
        if np.max(abs(remaining)) <= relative_tolerance:
            return trial_density
    raise ValueError(f"Grid cannot honor local density within {relative_tolerance:.1%}: "
                     f"worst residual {np.max(abs(remaining)):.2%}")


def candidate_grid(grid: IonosphereGrid, basis: TrackLogDensityBasis,
                   coefficients: np.ndarray, *,
                   local_density: LocalDensity | None = None,
                   maximum_log_correction: float = 0.5,
                   taper_depth_km: float = 200.0) -> IonosphereGrid:
    """Apply a bounded correction and optional in-situ constraint to a prior."""
    update = basis.correction(grid, coefficients)
    if np.max(abs(update)) > maximum_log_correction:
        raise ValueError("Trial exceeds the log-density trust region")
    density = np.asarray(grid.iono_en_grid, dtype=float) * np.exp(update)
    if local_density is not None:
        density = _anchor_upper_tail(grid, density, local_density, basis.frame,
                                     taper_depth_km=taper_depth_km,
                                     relative_tolerance=local_density.relative_uncertainty)
    if np.any(density <= 0) or not np.all(np.isfinite(density)):
        raise ValueError("Invalid candidate electron density")
    collision = grid.collision_freq
    if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
            and grid.neutral_species_cm3 is not None):
        collision = effective_collision_frequency(
            grid.electron_temp_k, grid.ion_temp_k, density * 1e6,
            grid.neutral_species_cm3)
    return replace(grid, iono_en_grid=density, iono_en_grid_5=density.copy(),
                   collision_freq=collision,
                   metadata={**grid.metadata, "retrieval": "general log-density field"})


def trust_region_step(baseline_score: float, plus_scores: np.ndarray,
                      minus_scores: np.ndarray, *,
                      perturbation: float, radius: float,
                      regularization: float = 1.0) -> np.ndarray:
    """Propose a bounded step from paired full-ray pilot scores.

    The score can include all accepted O/X returns and Doppler. A combined
    candidate still needs an independent full-ray evaluation before adoption.
    """
    plus = np.asarray(plus_scores, dtype=float)
    minus = np.asarray(minus_scores, dtype=float)
    if (plus.shape != minus.shape or not np.all(np.isfinite(plus))
            or not np.all(np.isfinite(minus)) or not np.isfinite(baseline_score)
            or perturbation <= 0 or radius <= 0 or regularization <= 0):
        raise ValueError("Invalid paired pilot scores or trust-region settings")
    gradient = (plus - minus) / (2.0 * perturbation)
    curvature = np.maximum((plus + minus - 2.0 * baseline_score)
                           / perturbation ** 2, 0.0)
    proposal = -gradient / (curvature + regularization)
    length = np.linalg.norm(proposal)
    return proposal * min(1.0, radius / length) if length > 0 else proposal
