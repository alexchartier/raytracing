from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np

from python_raytrace.grid import load_ionosphere_grid_netcdf
from reports.ionogram_metrics import Ionogram


ROOT = Path(__file__).resolve().parents[2]


class LatWavePassTests(unittest.TestCase):
    def test_pass_geometry_and_staged_density(self) -> None:
        manifest = json.loads((ROOT / "reports/data/lat_wave_pass_manifest.json").read_text())
        profiles = manifest["profiles"]
        latitude = np.array([row["latitude_deg"] for row in profiles])
        self.assertEqual(len(latitude), 20)
        self.assertAlmostEqual(latitude[-1] - latitude[0], 20.0)
        np.testing.assert_allclose(np.diff(latitude), 20.0 / 19)
        self.assertEqual(len({row["longitude_deg"] for row in profiles}), 1)

        with np.load(ROOT / manifest["density_grid"], allow_pickle=False) as truth:
            density = np.asarray(truth["electron_density_cm3"])
            altitudes = np.asarray(truth["altitudes_km"])
        forward = load_ionosphere_grid_netcdf(ROOT / manifest["forward_grid"])
        np.testing.assert_array_equal(forward.iono_en_grid, density)
        with np.load(ROOT / "reports/data/d_inverse_iri_truth_density.npz", allow_pickle=False) as baseline:
            original = np.asarray(baseline["electron_density_cm3"])
        near_peak = int(np.argmin(abs(altitudes - 280.0)))
        change = density[:, :, near_peak] / original[:, :, near_peak] - 1.0
        self.assertGreater(float(np.max(change)), 0.15)
        self.assertLess(float(np.min(change)), -0.15)

    def test_recovered_ionograms_keep_all_valid_reflected_returns(self) -> None:
        raw_dir = ROOT / "reports/data/lat_wave_pass_ionograms"
        recovered_dir = ROOT / "reports/data/lat_wave_pass_ionograms_recovered"
        recovered_total = 0
        for index in range(1, 21):
            raw = Ionogram.read(raw_dir / f"ionogram_{index:02d}.npz")
            final = Ionogram.read(recovered_dir / f"ionogram_{index:02d}.npz")
            with np.load(recovered_dir / f"ionogram_{index:02d}.npz", allow_pickle=False) as data:
                added = int(data["dense_recovered_return_count"])
                counts = np.asarray(data["count_array"], dtype=int)
            self.assertEqual(len(final.records), len(raw.records) + added)
            self.assertTrue(np.all(final.records[:, 2] >= 150.0))
            self.assertEqual(int(np.sum(counts)), len(final.records))
            recovered_total += added
        self.assertEqual(recovered_total, 43)


if __name__ == "__main__":
    unittest.main()
