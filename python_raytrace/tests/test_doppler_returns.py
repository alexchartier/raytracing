import numpy as np

from reports.doppler_returns import unique_returns


def test_duplicate_grouping_preserves_distinct_range_doppler_and_mode() -> None:
    records = np.array([
        [10, 1, 300.0],
        [10, 1, 300.3],  # Numerical duplicate of the first return.
        [10, 1, 300.4],  # Similar range, distinct Doppler branch.
        [10, 1, 303.0],  # Distinct range branch.
        [10, -1, 300.0],  # Distinct mode.
    ])
    doppler = np.array([2.0, 2.02, 3.0, 2.0, 2.0])
    unique = unique_returns(records, doppler)
    assert len(unique) == 4
    o_mode = unique[unique[:, 1] == 1]
    np.testing.assert_allclose(sorted(o_mode[:, 2]), [300.15, 300.4, 303.0])
    np.testing.assert_allclose(sorted(o_mode[:, 3]), [2.0, 2.01, 3.0])
