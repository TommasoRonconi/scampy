"""Tests for scampy.halo — mass functions, bias functions, density profiles.

All tests use the same cosmology and power-spectrum fixtures defined in
conftest.py.  The mass range [1e11, 4e12] M_⊙ is chosen where the mass
function is well-behaved for the simple P(k) ∝ k power law in simple_pk.

Note: the simple_pk fixture (P∝k) is unrealistic but sufficient to exercise
the code paths.  Qualitative tests (monotonicity, redshift dependence) are
therefore limited to the range and direction actually observed for this fixture.
All inputs are passed as arrays (not scalars) to avoid an indexing limitation
in the current implementation.
"""

import numpy as np
import pytest

from scampy.halo.mass_function import ShethTormen01, Tinker08
from scampy.halo.bias import ShethMoTormen01, Tinker10
from scampy.halo.density_profile import (
    concentration_shimizu03,
    density_profile_FT,
)


MM     = np.geomspace(1e11, 4e12, 8)   # safe range for P∝k; all values positive
MM_REF = np.array([1e11, 1e12, 4e12])  # 3-point regression reference


# ─────────────────────────────────────────────────────────────────────────────
# Mass functions
# ─────────────────────────────────────────────────────────────────────────────

class TestShethTormen01:

    def test_positive(self, simple_pk):
        assert np.all(ShethTormen01(MM, 0.0, simple_pk) > 0)

    def test_monotone_decreasing(self, simple_pk):
        """Halo abundance decreases with mass for this mass range."""
        n = ShethTormen01(MM, 0.0, simple_pk)
        assert np.all(np.diff(n) < 0)

    def test_output_shape(self, simple_pk):
        result = ShethTormen01(MM, 0.0, simple_pk)
        assert result.shape == MM.shape

    def test_regression_values(self, simple_pk):
        expected = np.array([1.26119363e-13, 2.63576217e-15, 1.98838629e-16])
        np.testing.assert_allclose(
            ShethTormen01(MM_REF, 0.0, simple_pk), expected, rtol=1e-4
        )


class TestTinker08:

    def test_positive(self, simple_pk):
        assert np.all(Tinker08(MM, 0.0, simple_pk) > 0)

    def test_monotone_decreasing(self, simple_pk):
        assert np.all(np.diff(Tinker08(MM, 0.0, simple_pk)) < 0)

    def test_output_shape(self, simple_pk):
        result = Tinker08(MM, 0.0, simple_pk)
        assert result.shape == MM.shape

    def test_tinker08_exceeds_st(self, simple_pk):
        """Tinker08 predicts higher abundance than ShethTormen01 throughout."""
        assert np.all(Tinker08(MM, 0.0, simple_pk) > ShethTormen01(MM, 0.0, simple_pk))

    def test_regression_values(self, simple_pk):
        expected = np.array([9.37283173e-13, 1.02809715e-14, 5.42080129e-16])
        np.testing.assert_allclose(
            Tinker08(MM_REF, 0.0, simple_pk), expected, rtol=1e-4
        )


# ─────────────────────────────────────────────────────────────────────────────
# Bias functions
# ─────────────────────────────────────────────────────────────────────────────

class TestShethMoTormen01:

    def test_positive(self, simple_pk):
        assert np.all(ShethMoTormen01(MM, 0.0, simple_pk) > 0)

    def test_output_shape(self, simple_pk):
        result = ShethMoTormen01(MM, 0.0, simple_pk)
        assert result.shape == MM.shape

    def test_finite(self, simple_pk):
        assert np.all(np.isfinite(ShethMoTormen01(MM, 0.0, simple_pk)))

    def test_regression_values(self, simple_pk):
        expected = np.array([0.9796463, 0.88893566, 0.77702419])
        np.testing.assert_allclose(
            ShethMoTormen01(MM_REF, 0.0, simple_pk), expected, rtol=1e-4
        )


class TestTinker10:

    def test_positive(self, simple_pk):
        assert np.all(Tinker10(MM, 0.0, simple_pk) > 0)

    def test_output_shape(self, simple_pk):
        result = Tinker10(MM, 0.0, simple_pk)
        assert result.shape == MM.shape

    def test_finite(self, simple_pk):
        assert np.all(np.isfinite(Tinker10(MM, 0.0, simple_pk)))

    def test_regression_values(self, simple_pk):
        expected = np.array([0.6470409, 0.60568823, 0.59160816])
        np.testing.assert_allclose(
            Tinker10(MM_REF, 0.0, simple_pk), expected, rtol=1e-4
        )


# ─────────────────────────────────────────────────────────────────────────────
# Concentration model and NFW profile FT
# ─────────────────────────────────────────────────────────────────────────────

class TestConcentrationShimizu03:

    def test_positive(self, simple_pk):
        cc = concentration_shimizu03(MM, 0.0, simple_pk)
        assert np.all(cc > 0)

    def test_decreases_with_mass(self, simple_pk):
        """More massive halos are less concentrated (c ∝ M^{-0.13})."""
        cc = concentration_shimizu03(MM, 0.0, simple_pk)
        assert np.all(np.diff(cc) < 0)

    def test_decreases_with_redshift(self, simple_pk):
        """c(M,z) = 8/(1+z) × f(M): higher z → lower concentration."""
        c0 = concentration_shimizu03(MM, 0.0, simple_pk)
        c1 = concentration_shimizu03(MM, 1.0, simple_pk)
        assert np.all(c0 > c1)

    def test_redshift_scaling(self, simple_pk):
        """c(M, z=1) = c(M, z=0)/2 because of the 1/(1+z) factor."""
        c0 = concentration_shimizu03(MM, 0.0, simple_pk)
        c1 = concentration_shimizu03(MM, 1.0, simple_pk)
        np.testing.assert_allclose(c1, c0 / 2.0, rtol=1e-10)

    def test_regression_value(self, simple_pk):
        """c(10¹² M⊙, z=0) = 8 × (1.0204×10⁻¹⁴)^{-0.13}."""
        m = np.array([1e12])
        expected = 8.0 * (1.0204e-14 * m[0]) ** (-0.13)
        c = concentration_shimizu03(m, 0.0, simple_pk)
        assert float(c) == pytest.approx(expected, rel=1e-6)


class TestDensityProfileFT:

    def test_finite_output(self, simple_pk):
        """NFW profile FT must be finite for all reasonable k."""
        m = 1e12
        kk = np.geomspace(0.01, 10.0, 20)
        result = density_profile_FT(kk, m, 0.0, simple_pk)
        assert np.all(np.isfinite(result))

    def test_approaches_one_at_k0(self, simple_pk):
        """Normalised NFW FT → 1 as k → 0."""
        m = 1e12
        kk = np.geomspace(1e-4, 0.01, 5)
        result = density_profile_FT(kk, m, 0.0, simple_pk)
        np.testing.assert_allclose(result, 1.0, atol=0.01)

    def test_decreases_with_k(self, simple_pk):
        """NFW profile FT is a decreasing function of k."""
        m = 1e12
        kk = np.geomspace(0.1, 10.0, 10)
        result = density_profile_FT(kk, m, 0.0, simple_pk)
        assert np.all(np.diff(result) <= 0)
