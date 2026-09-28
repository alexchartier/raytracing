"""Group numerical near-duplicate ionogram returns for Doppler scoring."""

from __future__ import annotations

import numpy as np


def unique_returns(records: np.ndarray, doppler: np.ndarray,
                   range_tolerance_km: float = 1.0,
                   doppler_tolerance_hz: float = 0.1) -> np.ndarray:
    """Return [frequency index, mode, group range, Doppler] for distinct paths.

    One accepted return is retained for each tight cluster. Distinct branches
    with different ranges or Doppler remain separate. The original accepted
    rays are never discarded from the saved ionograms.
    """
    records = np.asarray(records, dtype=float)
    doppler = np.asarray(doppler, dtype=float)
    if records.ndim != 2 or records.shape[1] < 3 or len(records) != len(doppler):
        raise ValueError("records and Doppler must be aligned")
    rows = []
    for mode in (1, -1):
        for frequency_index in range(81):
            indices = np.flatnonzero((records[:, 1] == mode)
                                    & (records[:, 0] == frequency_index))
            if not len(indices):
                continue
            ordered = indices[np.argsort(records[indices, 2])]
            clusters: list[list[int]] = []
            for index in ordered:
                if (clusters
                        and records[index, 2] - records[clusters[-1][0], 2]
                        <= range_tolerance_km
                        and abs(doppler[index] - np.median(doppler[clusters[-1]]))
                        <= doppler_tolerance_hz):
                    clusters[-1].append(int(index))
                else:
                    clusters.append([int(index)])
            for cluster in clusters:
                rows.append((float(frequency_index), float(mode),
                             float(np.median(records[cluster, 2])),
                             float(np.median(doppler[cluster]))))
    return np.asarray(rows, dtype=float).reshape(-1, 4)
