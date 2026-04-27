"""Tests for scampy.utilities.functions and scampy.utilities.constants."""

import numpy as np
import pytest

from scampy.utilities import constants
from scampy.utilities.functions import (
    FT_tophat,
    FT_tophat_D1,
    trap_int,
    linear_interpolation,
)


# ─────────────────────────────────────────────────────────────────────────────
# constants
# ─────────────────────────────────────────────────────────────────────────────

class TestConstants:

    def test_degrees_to_radians_value(self):
        assert constants.degrees_to_radians == pytest.approx(np.pi / 180.0, rel=1e-12)

    def test_degrees_to_radians_right_angle(self):
        assert 90.0 * constants.degrees_to_radians == pytest.approx(np.pi / 2.0, rel=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# FT_tophat
# ─────────────────────────────────────────────────────────────────────────────

class TestFTTophat:

    def test_small_kR_approaches_one(self):
        """W(kR) → 1 as kR → 0; kR=1e-3 avoids catastrophic cancellation."""
        val = FT_tophat(1e-3)
        assert val == pytest.approx(1.0, rel=1e-5)

    def test_positive_for_small_kR(self):
        kR = np.array([0.01, 0.1, 0.5])
        assert np.all(FT_tophat(kR) > 0)

    def test_range_kR_one(self):
        """At kR=1, W ≈ 0.9035 (analytical: 3(sin1-cos1))."""
        expected = 3.0 * (np.sin(1.0) - np.cos(1.0))
        assert FT_tophat(1.0) == pytest.approx(expected, rel=1e-6)

    def test_range_kR_pi(self):
        """At kR=π, W = 3(sin π - π cos π)/π^3 = 3π/π^3 = 3/π^2."""
        kR = np.pi
        expected = 3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR**3
        assert FT_tophat(kR) == pytest.approx(expected, rel=1e-10)

    def test_vectorised(self):
        kR = np.array([0.01, 0.1, 1.0, 3.0])
        result = FT_tophat(kR)
        assert result.shape == (4,)

    def test_regression_values(self):
        kR = np.array([0.01, 0.1, 1.0, 3.0])
        expected = np.array([0.99999, 0.99900036, 0.90350604, 0.3456775])
        np.testing.assert_allclose(FT_tophat(kR), expected, rtol=1e-5)


class TestFTTophatD1:

    def test_negative_for_small_kR(self):
        """Derivative W'(kR) is negative for small positive kR."""
        kR = np.array([0.01, 0.1, 1.0])
        assert np.all(FT_tophat_D1(kR) < 0)

    def test_regression_values(self):
        kR = np.array([0.01, 0.1, 1.0, 3.0])
        expected = np.array([-0.00199998, -0.01998572, -0.18610516, -0.2986375])
        np.testing.assert_allclose(FT_tophat_D1(kR), expected, rtol=1e-5)


# ─────────────────────────────────────────────────────────────────────────────
# trap_int
# ─────────────────────────────────────────────────────────────────────────────

class TestTrapInt:

    def test_constant_function(self):
        """∫₀¹ 3 dx = 3."""
        x = np.linspace(0.0, 1.0, 1001)
        assert trap_int(x, np.full_like(x, 3.0)) == pytest.approx(3.0, rel=1e-10)

    def test_linear_function(self):
        """∫₀¹ x dx = 0.5 (trapezoidal rule is exact for linear)."""
        x = np.linspace(0.0, 1.0, 1001)
        assert trap_int(x, x) == pytest.approx(0.5, rel=1e-10)

    def test_quadratic_function(self):
        """∫₀¹ x² dx = 1/3 (trapezoidal with many points)."""
        x = np.linspace(0.0, 1.0, 10001)
        assert trap_int(x, x**2) == pytest.approx(1.0 / 3.0, rel=1e-6)

    def test_exponential_function(self):
        """∫₀² e^x dx = e² - 1 ≈ 6.389."""
        x = np.linspace(0.0, 2.0, 10001)
        assert trap_int(x, np.exp(x)) == pytest.approx(np.e**2 - 1.0, rel=1e-6)

    def test_single_interval(self):
        """Two points: trapezoid is exact."""
        x = np.array([0.0, 1.0])
        y = np.array([2.0, 4.0])
        assert trap_int(x, y) == pytest.approx(3.0, rel=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# linear_interpolation
# ─────────────────────────────────────────────────────────────────────────────

class TestLinearInterpolation:

    def test_at_grid_nodes(self):
        """Querying at grid nodes must return exact values."""
        xg = np.array([0.0, 1.0, 2.0, 3.0])
        yg = np.array([0.0, 1.0, 4.0, 9.0])
        np.testing.assert_allclose(
            linear_interpolation(xg, xg, yg), yg, rtol=1e-12
        )

    def test_at_midpoints_linear_function(self):
        """Piecewise linear: midpoints of y=x are exact."""
        xg = np.array([0.0, 1.0, 2.0, 3.0])
        yg = xg.copy()
        xq = np.array([0.5, 1.5, 2.5])
        np.testing.assert_allclose(
            linear_interpolation(xq, xg, yg), xq, rtol=1e-12
        )

    def test_regression_quadratic(self):
        """Midpoints of y=x² are not exact but have specific values."""
        xg = np.array([0.0, 1.0, 2.0, 3.0])
        yg = xg**2
        xq = np.array([0.5, 1.5, 2.5])
        expected = np.array([0.5, 2.5, 6.5])   # linear interp of parabola
        np.testing.assert_allclose(
            linear_interpolation(xq, xg, yg), expected, rtol=1e-12
        )

    def test_output_shape(self):
        xg = np.linspace(0, 10, 20)
        yg = np.sin(xg)
        xq = np.linspace(1, 9, 50)
        result = linear_interpolation(xq, xg, yg)
        assert result.shape == (50,)
