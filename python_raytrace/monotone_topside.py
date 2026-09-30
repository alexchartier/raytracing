"""Spatially varying, monotone topside profiles constrained by in-situ density.

The shape at each along-track position is specified by two positive log-density
fractions between its F2 peak and the spacecraft. PCHIP joins those fractions
without an overshoot. An in-situ measurement at every kilometre is checked
against the resulting forward-grid interpolation.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator

from .general_field_inverse import LocalDensity
from .grid import IonosphereGrid


def monotone_topside_grid(
    grid: IonosphereGrid,
    local: LocalDensity,
    sounding_latitudes_deg: np.ndarray,
    fraction_at_quarter: np.ndarray,
    fraction_at_three_fifths: np.ndarray,
    peak_height_shift_km: np.ndarray,
) -> IonosphereGrid:
    """Apply latitude-varying height and topside shape, fixing every 800 km sample.

    Fractions are of the log-density span from local F2 peak to measured
    spacecraft density: 1 at the peak and 0 at the spacecraft. Hence
    ``1 > fraction_at_quarter > fraction_at_three_fifths > 0`` guarantees a
    falling topside. A shift moves only the peak and bottomside; the topside
    is newly constructed from the shifted peak to spacecraft altitude.
    """
    latitudes = np.asarray(sounding_latitudes_deg, dtype=float)
    q1 = np.asarray(fraction_at_quarter, dtype=float)
    q2 = np.asarray(fraction_at_three_fifths, dtype=float)
    shift = np.asarray(peak_height_shift_km, dtype=float)
    if (latitudes.ndim != 1 or len(latitudes) < 2
            or not np.all(np.diff(latitudes) > 0)
            or any(value.shape != latitudes.shape for value in (q1, q2, shift))
            or not np.all(np.isfinite(np.stack((q1, q2, shift))))
            or np.any(q2 <= 0) or np.any(q1 >= 1) or np.any(q1 <= q2)
            or np.max(abs(shift)) > 60):
        raise ValueError("Invalid spatial topside parameters")
    points = np.asarray(local.lat_lon_deg, dtype=float)
    if (np.ptp(points[:, 1]) > 1e-8 or
            np.max(np.diff(points[:, 0])) * 111.2 > 1.01 or
            np.min(points[:, 0]) > latitudes[0] or
            np.max(points[:, 0]) < latitudes[-1]):
        raise ValueError("Expected a dense, ordered along-track in-situ pass")
    alt = np.asarray(grid.altitudes_km, dtype=float)
    lon = np.asarray(grid.longitudes_deg, dtype=float)
    if (local.altitude_km not in alt or
            not lon[0] <= points[0, 1] <= lon[-1]):
        raise ValueError("The spacecraft altitude must be a grid node and the track longitude must be inside the grid")
    sc_index = int(np.where(alt == local.altitude_km)[0][0])
    original = np.asarray(grid.iono_en_grid, dtype=float)
    base_at_sc = original[:, :, sc_index]
    track_lat = points[:, 0]
    track_density = np.asarray(local.electron_density_cm3, dtype=float)
    if not np.all(np.diff(track_lat) > 0):
        raise ValueError("The in-situ track must be sorted by latitude")

    # At the measured longitude, recover grid-node densities from the 1 km
    # observations. The saved SAMI observation is piecewise linear on these
    # same latitude nodes; interpolation then honors every 1 km sample.
    target_at_nodes = np.interp(np.clip(grid.latitudes_deg, track_lat[0], track_lat[-1]),
                                track_lat, track_density)
    # If a spacecraft track starts or ends between latitude grid nodes, its
    # first/last observed interval needs the enclosing node as well. Extend
    # the measured local slope to that one node only, preserving the linear
    # interpolation of every measurement inside the track.
    below = np.flatnonzero(grid.latitudes_deg < track_lat[0])
    above = np.flatnonzero(grid.latitudes_deg > track_lat[-1])
    if len(below):
        edge = below[-1]
        slope = (track_density[1] - track_density[0]) / (track_lat[1] - track_lat[0])
        target_at_nodes[edge] = track_density[0] + slope * (grid.latitudes_deg[edge] - track_lat[0])
    if len(above):
        edge = above[0]
        slope = (track_density[-1] - track_density[-2]) / (track_lat[-1] - track_lat[-2])
        target_at_nodes[edge] = track_density[-1] + slope * (grid.latitudes_deg[edge] - track_lat[-1])
    if np.any(target_at_nodes <= 0):
        raise ValueError("Track-edge extrapolation produced nonpositive spacecraft density")
    base_at_track = np.array([np.interp(points[0, 1], lon, row)
                              for row in base_at_sc])
    correction = target_at_nodes / base_at_track
    target_at_sc = base_at_sc * correction[:, None]
    sample_lat = np.clip(np.asarray(grid.latitudes_deg, dtype=float),
                         latitudes[0], latitudes[-1])
    latitude_q1 = PchipInterpolator(latitudes, q1)(sample_lat)
    latitude_q2 = PchipInterpolator(latitudes, q2)(sample_lat)
    latitude_shift = PchipInterpolator(latitudes, shift)(sample_lat)
    revised = np.empty_like(original)
    for i in range(len(grid.latitudes_deg)):
        for j in range(len(lon)):
            old = original[i, j]
            old_peak = int(np.argmax(old))
            peak_height = alt[old_peak] + latitude_shift[i]
            peak_density = float(old[old_peak])
            spacecraft_density = float(target_at_sc[i, j])
            if (peak_height <= alt[0] or peak_height >= local.altitude_km - 100
                    or spacecraft_density <= 0 or spacecraft_density >= peak_density):
                raise ValueError("Peak and in-situ anchors are incompatible")
            low = alt < peak_height
            revised[i, j, low] = np.interp(alt[low] - latitude_shift[i], alt, old)
            span = np.log(peak_density / spacecraft_density)
            knots = peak_height + np.array([0.0, .25, .60, 1.0]) * (
                local.altitude_km - peak_height)
            log_values = np.log(spacecraft_density) + span * np.array(
                [1.0, latitude_q1[i], latitude_q2[i], 0.0])
            middle = (alt >= peak_height) & (alt <= local.altitude_km)
            revised[i, j, middle] = np.exp(PchipInterpolator(knots, log_values)(alt[middle]))
            upper = alt > local.altitude_km
            tail_slope = -span * latitude_q2[i] / (
                .4 * (local.altitude_km - peak_height))
            revised[i, j, upper] = spacecraft_density * np.exp(
                tail_slope * (alt[upper] - local.altitude_km))
            peak = int(np.argmax(revised[i, j]))
            if np.any(np.diff(revised[i, j, peak:]) > 1e-8 * peak_density):
                raise ValueError("Generated topside is not monotone")
    result = replace(grid, iono_en_grid=revised, iono_en_grid_5=revised.copy(),
                     metadata={**grid.metadata,
                               "retrieval": "spatially varying monotone topside"})
    sampled = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, alt), revised,
        bounds_error=True)(np.column_stack((points,
                                            np.full(len(points), local.altitude_km))))
    error = np.max(abs(sampled / track_density - 1))
    if error > local.relative_uncertainty:
        raise ValueError(f"Dense in-situ constraint missed by {100 * error:.3f}%")
    return result
