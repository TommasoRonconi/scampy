"""Tests for scampy.boxes — coordinate transforms and box utilities.

All reference values computed from the heritage C++ implementation.
Pure-math functions; no fixtures required.
"""

import numpy as np
import pytest

from scampy.boxes import (
    shift_periodic_box,
    centre_box,
    cartesian_to_polar,
    polar_to_cartesian,
    cartesian_to_equatorial,
    equatorial_to_cartesian,
    angular_to_euclidean_dist,
    random_projection,
)

# Fixed 3-point catalogue (Nobj, Ndim) used across tests.
_BOX = np.array([[1.0, 2.0, 3.0],
                 [99.0, 99.0, 99.0],
                 [50.0, 50.0, 50.0]])
_LBOX = 100.0


# ─────────────────────────────────────────────────────────────────────────────
# shift_periodic_box
# ─────────────────────────────────────────────────────────────────────────────

class TestShiftPeriodicBox:

    def test_scalar_shift_regression(self):
        result = shift_periodic_box(_BOX, 10.0, _LBOX)
        expected = np.array([[11., 12., 13.], [9., 9., 9.], [60., 60., 60.]])
        np.testing.assert_array_equal(result, expected)

    def test_vector_shift_regression(self):
        result = shift_periodic_box(_BOX, [5.0, -10.0, 20.0], _LBOX)
        expected = np.array([[6., 92., 23.], [4., 89., 19.], [55., 40., 70.]])
        np.testing.assert_array_equal(result, expected)

    def test_output_in_bounds(self):
        """All coordinates must stay in [0, Lbox) after shift."""
        rng = np.random.default_rng(1)
        box = rng.uniform(0, _LBOX, (100, 3))
        result = shift_periodic_box(box, 37.5, _LBOX)
        assert np.all(result >= 0.0)
        assert np.all(result <= _LBOX)

    def test_does_not_mutate_input(self):
        box = _BOX.copy()
        _ = shift_periodic_box(box, 10.0, _LBOX)
        np.testing.assert_array_equal(box, _BOX)

    def test_zero_shift_identity(self):
        result = shift_periodic_box(_BOX, 0.0, _LBOX)
        np.testing.assert_array_equal(result, _BOX)

    def test_wrong_shift_length_raises(self):
        with pytest.raises(AttributeError):
            shift_periodic_box(_BOX, [1.0, 2.0], _LBOX)   # Ndim mismatch

    def test_full_box_shift_is_identity(self):
        """Shifting by exactly Lbox wraps all coordinates back to themselves."""
        result = shift_periodic_box(_BOX, _LBOX, _LBOX)
        np.testing.assert_allclose(result, _BOX, atol=1e-10)


# ─────────────────────────────────────────────────────────────────────────────
# centre_box
# ─────────────────────────────────────────────────────────────────────────────

class TestCentreBox:

    def test_regression(self):
        result = centre_box(_BOX, _LBOX)
        expected = np.array([[-49., -48., -47.], [49., 49., 49.], [0., 0., 0.]])
        np.testing.assert_array_equal(result, expected)

    def test_box_centre_maps_to_zero(self):
        """A point at Lbox/2 in every direction should map to (0,0,0)."""
        centre = np.array([[_LBOX / 2, _LBOX / 2, _LBOX / 2]])
        result = centre_box(centre, _LBOX)
        np.testing.assert_array_equal(result, np.zeros((1, 3)))

    def test_does_not_mutate_input(self):
        box = _BOX.copy()
        _ = centre_box(box, _LBOX)
        np.testing.assert_array_equal(box, _BOX)


# ─────────────────────────────────────────────────────────────────────────────
# cartesian_to_polar / polar_to_cartesian
# ─────────────────────────────────────────────────────────────────────────────

# Three test points: (1,0,0), (0,1,0), (-1,0,1)
_COORDS = np.array([[1.0, 0.0, -1.0],
                    [0.0, 1.0,  0.0],
                    [0.0, 0.0,  1.0]])


class TestCartesianPolar:

    def test_rho_regression(self):
        rho, _, _ = cartesian_to_polar(_COORDS)
        expected = np.array([1.0, 1.0, np.sqrt(2)])
        np.testing.assert_allclose(rho, expected, rtol=1e-10)

    def test_theta_regression(self):
        _, theta, _ = cartesian_to_polar(_COORDS)
        expected = np.array([np.pi / 2, np.pi / 2, np.pi / 4])
        np.testing.assert_allclose(theta, expected, rtol=1e-10)

    def test_phi_regression(self):
        _, _, phi = cartesian_to_polar(_COORDS)
        expected = np.array([0.0, np.pi / 2, np.pi])
        np.testing.assert_allclose(phi, expected, rtol=1e-10)

    def test_round_trip(self):
        """polar_to_cartesian(cartesian_to_polar(x)) = x."""
        rho, theta, phi = cartesian_to_polar(_COORDS)
        x2, y2, z2 = polar_to_cartesian((rho, theta, phi))
        np.testing.assert_allclose(x2, _COORDS[0], atol=1e-10)
        np.testing.assert_allclose(y2, _COORDS[1], atol=1e-10)
        np.testing.assert_allclose(z2, _COORDS[2], atol=1e-10)

    def test_with_centre(self):
        """Non-zero centre shifts the origin correctly."""
        centre = (1.0, 0.0, 0.0)
        rho, _, _ = cartesian_to_polar(_COORDS, centre=centre)
        # point (1,0,0) relative to centre (1,0,0) has rho=0
        assert rho[0] == pytest.approx(0.0, abs=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# cartesian_to_equatorial / equatorial_to_cartesian
# ─────────────────────────────────────────────────────────────────────────────

class TestCartesianEquatorial:

    def test_d_regression(self):
        dd, _, _ = cartesian_to_equatorial(_COORDS)
        expected = np.array([1.0, 1.0, np.sqrt(2)])
        np.testing.assert_allclose(dd, expected, rtol=1e-10)

    def test_ra_regression(self):
        _, ra, _ = cartesian_to_equatorial(_COORDS)
        expected = np.array([0.0, np.pi / 2, np.pi])
        np.testing.assert_allclose(ra, expected, rtol=1e-10)

    def test_dec_regression(self):
        _, _, dec = cartesian_to_equatorial(_COORDS)
        expected = np.array([0.0, 0.0, np.pi / 4])
        np.testing.assert_allclose(dec, expected, rtol=1e-10)

    def test_round_trip(self):
        """equatorial_to_cartesian(cartesian_to_equatorial(x)) = x."""
        dd, ra, dec = cartesian_to_equatorial(_COORDS)
        x3, y3, z3 = equatorial_to_cartesian((dd, ra, dec))
        np.testing.assert_allclose(x3, _COORDS[0], atol=1e-10)
        np.testing.assert_allclose(y3, _COORDS[1], atol=1e-10)
        np.testing.assert_allclose(z3, _COORDS[2], atol=1e-10)


# ─────────────────────────────────────────────────────────────────────────────
# angular_to_euclidean_dist
# ─────────────────────────────────────────────────────────────────────────────

class TestAngularToEuclideanDist:

    def test_zero_angle_gives_zero(self):
        assert angular_to_euclidean_dist(0.0) == 0.0

    def test_pi_over_2_gives_sqrt2_r(self):
        """Angular distance π/2 on unit sphere gives chord √2."""
        assert angular_to_euclidean_dist(np.pi / 2) == pytest.approx(np.sqrt(2), rel=1e-10)

    def test_regression_array(self):
        d = np.array([0.0, np.pi / 6, np.pi / 3, np.pi / 2])
        expected = np.array([0.0, 0.51763809, 1.0, np.sqrt(2)])
        np.testing.assert_allclose(angular_to_euclidean_dist(d), expected, rtol=1e-6)

    def test_radius_scaling(self):
        """Result scales linearly with sphere radius r."""
        d = np.array([np.pi / 6, np.pi / 3, np.pi / 2])
        r1 = angular_to_euclidean_dist(d, r=1.0)
        r2 = angular_to_euclidean_dist(d, r=2.0)
        np.testing.assert_allclose(r2, 2.0 * r1, rtol=1e-12)

    def test_regression_r2(self):
        d = np.array([0.0, np.pi / 6, np.pi / 3, np.pi / 2])
        expected = np.array([0.0, 1.03527618, 2.0, 2.82842712])
        np.testing.assert_allclose(angular_to_euclidean_dist(d, r=2.0), expected, rtol=1e-6)


# ─────────────────────────────────────────────────────────────────────────────
# random_projection
# ─────────────────────────────────────────────────────────────────────────────

class TestRandomProjection:

    def test_output_shapes(self):
        RA, Dec = random_projection((0.0, 1.0), (-0.5, 0.5), size=20, rng=0)
        assert RA.shape == (20,)
        assert Dec.shape == (20,)

    def test_ra_in_limits(self):
        RA_lim = (0.1, 0.9)
        RA, _ = random_projection(RA_lim, (-0.4, 0.4), size=200, rng=1)
        assert np.all(RA >= RA_lim[0]) and np.all(RA <= RA_lim[1])

    def test_dec_in_limits(self):
        Dec_lim = (-0.3, 0.3)
        _, Dec = random_projection((0.0, 1.0), Dec_lim, size=200, rng=2)
        assert np.all(Dec >= Dec_lim[0]) and np.all(Dec <= Dec_lim[1])

    def test_regression_seed42(self):
        RA, Dec = random_projection((0.0, 1.0), (-0.5, 0.5), size=5, rng=42)
        expected_RA  = np.array([0.97562235, 0.7611397, 0.78606431, 0.12811363, 0.45038594])
        expected_Dec = np.array([0.26580186, -0.05864008, 0.35100532, 0.19039477, -0.39967990])
        np.testing.assert_allclose(RA,  expected_RA,  rtol=1e-6)
        np.testing.assert_allclose(Dec, expected_Dec, rtol=1e-6)

    def test_size_one_returns_scalar_like(self):
        """size=1 returns arrays of length 1."""
        RA, Dec = random_projection((0.0, 1.0), (-0.5, 0.5), size=1, rng=3)
        assert np.asarray(RA).size == 1
        assert np.asarray(Dec).size == 1
