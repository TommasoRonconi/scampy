"""Tests for scampy.lightcone — angular sizes, projections, lightcone limits.

Reference values computed from the heritage C++ implementation.
"""

import numpy as np
import pytest

from scampy.lightcone import (
    sequential_lightcones_limits,
    ang_size,
    size_ang,
    max_fov,
)
from scampy.utilities.interpolation import lin_interp


Z_ARR = np.array([0.1, 0.5, 1.0, 2.0])


@pytest.fixture(scope="module")
def d2z(default_cosmo):
    """Distance-to-redshift interpolator built from default_cosmo.dC."""
    zz = np.linspace(0.0, 20.0, 10000)
    dc = default_cosmo.dC(zz)
    return lin_interp(dc, zz)


# ─────────────────────────────────────────────────────────────────────────────
# ang_size  — comoving size → angular size
# ─────────────────────────────────────────────────────────────────────────────

class TestAngSize:

    def test_positive(self, default_cosmo):
        assert ang_size(1.0, 0.5, default_cosmo) > 0

    def test_scalar_regression(self, default_cosmo):
        """ang_size(1 Mpc/h, z=0.5) matches heritage value."""
        assert ang_size(1.0, 0.5, default_cosmo) == pytest.approx(5.2949e-04, rel=1e-3)

    def test_array_regression(self, default_cosmo):
        expected = np.array([0.00238975, 0.00052949, 0.00030268, 0.00019305])
        result = ang_size(1.0, Z_ARR, default_cosmo)
        np.testing.assert_allclose(result, expected, rtol=1e-3)

    def test_decreasing_with_z(self, default_cosmo):
        """For fixed physical size, angular size decreases with distance."""
        angles = ang_size(1.0, Z_ARR, default_cosmo)
        assert np.all(np.diff(angles) < 0)

    def test_output_shape(self, default_cosmo):
        result = ang_size(1.0, Z_ARR, default_cosmo)
        assert result.shape == Z_ARR.shape

    def test_proportional_to_size(self, default_cosmo):
        """Doubling lC at fixed z must double the angular size."""
        a1 = ang_size(1.0, 0.5, default_cosmo)
        a2 = ang_size(2.0, 0.5, default_cosmo)
        assert a2 == pytest.approx(2.0 * a1, rel=1e-10)


# ─────────────────────────────────────────────────────────────────────────────
# size_ang  — angular size → comoving size
# ─────────────────────────────────────────────────────────────────────────────

class TestSizeAng:

    def test_positive(self, default_cosmo):
        assert size_ang(1e-3, 0.5, default_cosmo) > 0

    def test_scalar_regression(self, default_cosmo):
        """size_ang(1e-3 rad, z=0.5) matches heritage value."""
        assert size_ang(1e-3, 0.5, default_cosmo) == pytest.approx(1.8886, rel=1e-3)

    def test_array_regression(self, default_cosmo):
        expected = np.array([0.41845449, 1.88862565, 3.30383249, 5.17988345])
        result = size_ang(1e-3, Z_ARR, default_cosmo)
        np.testing.assert_allclose(result, expected, rtol=1e-3)

    def test_increasing_with_z(self, default_cosmo):
        """Fixed angular size corresponds to larger physical size at higher z."""
        sizes = size_ang(1e-3, Z_ARR, default_cosmo)
        assert np.all(np.diff(sizes) > 0)

    def test_output_shape(self, default_cosmo):
        result = size_ang(1e-3, Z_ARR, default_cosmo)
        assert result.shape == Z_ARR.shape

    def test_round_trip_scalar(self, default_cosmo):
        """size_ang(ang_size(lC, z)) = lC."""
        for lC in [1.0, 5.0, 10.0]:
            th = ang_size(lC, 0.5, default_cosmo)
            lC_back = size_ang(th, 0.5, default_cosmo)
            assert lC_back == pytest.approx(lC, rel=1e-10)

    def test_round_trip_array(self, default_cosmo):
        """Vectorised round-trip: size_ang(ang_size(lC, z_arr)) = lC."""
        lC_arr = np.array([1.0, 5.0, 10.0])
        th_arr = ang_size(lC_arr, Z_ARR[:3], default_cosmo)
        lC_back = size_ang(th_arr, Z_ARR[:3], default_cosmo)
        np.testing.assert_allclose(lC_back, lC_arr, rtol=1e-10)


# ─────────────────────────────────────────────────────────────────────────────
# max_fov  — solid angle of the field of view
# ─────────────────────────────────────────────────────────────────────────────

class TestMaxFov:

    def test_regression_square_fov(self):
        """A square FOV with dxmax=dymax=2*rad fills a full circle: π sr."""
        result = max_fov(rad=250.0, dxmax=500.0, dymax=500.0)
        assert result == pytest.approx(np.pi, rel=1e-6)

    def test_regression_rectangular_fov(self):
        result = max_fov(rad=250.0, dxmax=300.0, dymax=400.0)
        assert result == pytest.approx(1.2870022175865685, rel=1e-5)

    def test_positive(self):
        assert max_fov(100.0, 200.0, 200.0) > 0

    def test_smaller_box_gives_smaller_fov(self):
        """A smaller box dimension must give a smaller or equal field of view."""
        fov_large = max_fov(250.0, 500.0, 500.0)
        fov_small = max_fov(250.0, 300.0, 300.0)
        assert fov_large >= fov_small

    def test_symmetric_in_dxmax_dymax(self):
        """Swapping dxmax and dymax must not change the result for equal values."""
        fov1 = max_fov(200.0, 300.0, 400.0)
        fov2 = max_fov(200.0, 400.0, 300.0)
        assert fov1 == pytest.approx(fov2, rel=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# sequential_lightcones_limits
# ─────────────────────────────────────────────────────────────────────────────

class TestSequentialLightconesLimits:

    def test_returns_two_arrays(self, d2z):
        result = sequential_lightcones_limits(zmax=1.0, Lbox=500.0, d2z=d2z)
        assert len(result) == 2
        limits, shifts = result
        assert isinstance(np.asarray(limits), np.ndarray)
        assert isinstance(np.asarray(shifts), np.ndarray)

    def test_shifts_non_negative(self, d2z):
        _, shifts = sequential_lightcones_limits(zmax=1.0, Lbox=500.0, d2z=d2z)
        assert np.all(np.asarray(shifts) >= 0.0)

    def test_first_shift_is_half_box(self, d2z):
        """With centre=True the first box is centred at Lbox/2."""
        _, shifts = sequential_lightcones_limits(
            zmax=1.0, Lbox=500.0, d2z=d2z, centre=True
        )
        assert np.asarray(shifts)[0] == pytest.approx(250.0)

    def test_last_limit_exceeds_zmax(self, d2z):
        """The loop exits only when the last redshift exceeds zmax."""
        limits, _ = sequential_lightcones_limits(zmax=1.0, Lbox=500.0, d2z=d2z)
        assert limits[-1] > 1.0

    def test_regression(self, d2z):
        limits, shifts = sequential_lightcones_limits(
            zmax=1.0, Lbox=500.0, d2z=d2z
        )
        expected_limits = np.array([0.05917114, 0.18288157, 0.31520124, 0.45816433,
                                    0.61420953, 0.78629523, 0.97801941, 1.19382816])
        expected_shifts = np.array([250., 750., 1250., 1750., 2250., 2750., 3250., 3750.])
        np.testing.assert_allclose(np.asarray(limits), expected_limits, rtol=1e-4)
        np.testing.assert_allclose(np.asarray(shifts), expected_shifts, rtol=1e-6)
