from __future__ import annotations

import datetime as dt
from dataclasses import dataclass

import numpy as np
import pymsis


@dataclass
class MsisAtmosphere:
    neutral_species_cm3: np.ndarray
    neutral_temperature_k: np.ndarray


_ITIKAWA_TEMPERATURE_AXIS_K = np.array(
    [100.0, 200.0, 300.0, 400.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0],
    dtype=float,
)

_ITIKAWA_COLLISION_DATA = np.array(
    [
        [0.255, 0.491, 0.723, 0.948, 1.17, 2.11, 2.86, 3.49, 4.08, 4.72, 5.44, 6.23, 7.08, 7.96],  # N2
        [0.123, 0.237, 0.335, 0.425, 0.512, 0.956, 1.41, 1.87, 2.32, 2.72, 3.08, 3.40, 3.68, 3.92],  # O2
        [0.589, 0.528, 0.566, 0.673, 0.859, 2.73, 4.88, 6.46, 7.49, 8.13, 8.54, 8.79, 8.95, 9.06],  # NO
        [158.0, 107.0, 84.9, 71.5, 62.3, 39.3, 29.3, 23.7, 20.1, 17.6, 15.8, 14.5, 13.4, 12.5],  # H2O
        [10.1, 10.1, 10.1, 9.96, 9.76, 8.04, 6.45, 5.34, 4.61, 4.14, 3.86, 3.71, 3.68, 3.74],  # CO2
        [1.51, 0.987, 0.719, 0.568, 0.480, 0.396, 0.486, 0.636, 0.822, 1.04, 1.29, 1.58, 1.91, 2.26],  # CH4
        [3.44, 4.85, 5.89, 6.72, 7.42, 9.81, 11.3, 12.2, 12.9, 13.4, 13.8, 14.1, 14.3, 14.5],  # H
        [0.448, 0.654, 0.820, 0.963, 1.09, 1.62, 2.05, 2.41, 2.74, 3.03, 3.30, 3.55, 3.78, 4.00],  # He
        [0.081, 0.138, 0.192, 0.245, 0.297, 0.551, 0.800, 1.04, 1.28, 1.51, 1.73, 1.94, 2.14, 2.34],  # O
        [0.341, 0.293, 0.243, 0.203, 0.173, 0.110, 0.126, 0.185, 0.272, 0.379, 0.504, 0.643, 0.795, 0.960],  # Ar
    ],
    dtype=float,
) * 1e-8

_PHARLAP_NEUTRAL_SPECIES_INDEX = {
    "N2": 0,
    "O2": 1,
    "H": 6,
    "HE": 7,
    "O": 8,
    "AR": 9,
}


def build_msis_atmosphere(
    when: dt.datetime,
    latitudes_deg: np.ndarray,
    longitudes_deg: np.ndarray,
    altitudes_km: np.ndarray,
    *,
    f107: float,
    f107a: float | None = None,
    ap_vector: np.ndarray | None = None,
    ap_daily: float | None = None,
    version: float | str = 2.1,
) -> MsisAtmosphere:
    if f107a is None:
        f107a = f107
    if ap_vector is None:
        ap_fill = 4.0 if ap_daily is None else float(ap_daily)
        ap_vector = np.full((7,), ap_fill, dtype=float)
    else:
        ap_vector = np.asarray(ap_vector, dtype=float)
        if ap_vector.shape != (7,):
            raise ValueError("ap_vector must have shape (7,)")
        ap_vector = ap_vector.copy()
        if ap_daily is not None:
            ap_vector[0] = float(ap_daily)

    output = pymsis.calculate(
        np.array([when], dtype=object),
        np.asarray(longitudes_deg, dtype=float),
        np.asarray(latitudes_deg, dtype=float),
        np.asarray(altitudes_km, dtype=float),
        f107s=np.array([float(f107)], dtype=float),
        f107as=np.array([float(f107a)], dtype=float),
        aps=np.asarray([ap_vector], dtype=float),
        version=version,
    )
    output_lat_lon_alt = np.transpose(output[0], (1, 0, 2, 3))

    neutral_species_cm3 = np.zeros((7, len(latitudes_deg), len(longitudes_deg), len(altitudes_km)), dtype=float)
    neutral_species_cm3[0] = output_lat_lon_alt[..., pymsis.Variable.HE] / 1e6
    neutral_species_cm3[1] = output_lat_lon_alt[..., pymsis.Variable.O] / 1e6
    neutral_species_cm3[2] = output_lat_lon_alt[..., pymsis.Variable.N2] / 1e6
    neutral_species_cm3[3] = output_lat_lon_alt[..., pymsis.Variable.O2] / 1e6
    neutral_species_cm3[4] = output_lat_lon_alt[..., pymsis.Variable.AR] / 1e6
    neutral_species_cm3[6] = output_lat_lon_alt[..., pymsis.Variable.H] / 1e6

    neutral_temperature_k = np.asarray(output_lat_lon_alt[..., pymsis.Variable.TEMPERATURE], dtype=float)
    return MsisAtmosphere(
        neutral_species_cm3=neutral_species_cm3,
        neutral_temperature_k=neutral_temperature_k,
    )


def electron_neutral_collision_frequency(
    electron_temperature_k: np.ndarray,
    neutral_number_density_cm3: np.ndarray,
    *,
    species_index: int,
) -> np.ndarray:
    temperature = np.asarray(electron_temperature_k, dtype=float)
    neutral_density = np.asarray(neutral_number_density_cm3, dtype=float)
    coeff = np.interp(
        np.clip(temperature, _ITIKAWA_TEMPERATURE_AXIS_K[0], _ITIKAWA_TEMPERATURE_AXIS_K[-1]),
        _ITIKAWA_TEMPERATURE_AXIS_K,
        _ITIKAWA_COLLISION_DATA[species_index],
        left=_ITIKAWA_COLLISION_DATA[species_index, 0],
        right=_ITIKAWA_COLLISION_DATA[species_index, -1],
    )
    return coeff * np.maximum(neutral_density, 0.0)


def electron_ion_collision_frequency(
    electron_temperature_k: np.ndarray,
    ion_temperature_k: np.ndarray,
    electron_density_m3: np.ndarray,
) -> np.ndarray:
    te = np.maximum(np.asarray(electron_temperature_k, dtype=float), 1.0)
    ti = np.maximum(np.asarray(ion_temperature_k, dtype=float), 1.0)
    ne = np.maximum(np.asarray(electron_density_m3, dtype=float), 0.0)

    out = np.zeros_like(ne, dtype=float)
    valid = ne > 0.0
    if not np.any(valid):
        return out

    te_valid = te[valid]
    ti_valid = ti[valid]
    ne_valid = ne[valid]

    ki_sq = 2.09985255e-4 * ne_valid / ti_valid
    ke_sq = 2.09985255e-4 * ne_valid / te_valid
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_coulomb = (
            13.484870477617616
            + np.log(te_valid)
            - 0.5 * np.log(ke_sq)
            - ((ke_sq + ki_sq) / ki_sq) * (0.5 * np.log((ki_sq + ke_sq) / ke_sq))
        )
        out[valid] = 3.63315e-6 * ne_valid * np.power(te_valid, -1.5) * ln_coulomb

    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def effective_collision_frequency(
    electron_temperature_k: np.ndarray,
    ion_temperature_k: np.ndarray,
    electron_density_m3: np.ndarray,
    neutral_species_cm3: np.ndarray,
) -> np.ndarray:
    neutrals = np.asarray(neutral_species_cm3, dtype=float)
    if neutrals.shape[0] != 7:
        raise ValueError("neutral_species_cm3 must use the PHaRLAP species ordering with leading dimension 7")

    total = (
        electron_neutral_collision_frequency(electron_temperature_k, neutrals[2], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["N2"])
        + electron_neutral_collision_frequency(electron_temperature_k, neutrals[3], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["O2"])
        + electron_neutral_collision_frequency(electron_temperature_k, neutrals[1], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["O"])
        + electron_neutral_collision_frequency(electron_temperature_k, neutrals[6], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["H"])
        + electron_neutral_collision_frequency(electron_temperature_k, neutrals[0], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["HE"])
        + electron_neutral_collision_frequency(electron_temperature_k, neutrals[4], species_index=_PHARLAP_NEUTRAL_SPECIES_INDEX["AR"])
        + electron_ion_collision_frequency(electron_temperature_k, ion_temperature_k, electron_density_m3)
    )
    return np.nan_to_num(total, nan=0.0, posinf=0.0, neginf=0.0)
