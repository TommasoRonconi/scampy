"""Tests for scampy.power_spectrum.power_spectrum.

All tests use the simple_pk fixture (P∝k) defined in conftest.py.
Reference values were computed from the heritage C++ implementation.
"""

import numpy as np
import pytest

from scampy.power_spectrum import power_spectrum


MM = np.array([1e11, 1e12, 4e12])
KK = np.array([0.01, 0.1, 1.0, 10.0])
RR = np.array([1.0, 5.0, 20.0])


# ─────────────────────────────────────────────────────────────────────────────
# Construction
# ─────────────────────────────────────────────────────────────────────────────

class TestConstruction:

    def test_default_cosmo_built(self):
        """Passing cosmo=None should construct a default cosmology."""
        kh = np.geomspace(1e-4, 1e2, 50)
        ps = power_spectrum(kh, kh)   # cosmo=None
        assert ps.cosmo is not None

    def test_wrong_cosmo_raises(self):
        kh = np.geomspace(1e-4, 1e2, 50)
        with pytest.raises(TypeError):
            power_spectrum(kh, kh, cosmo="not a model")

    def test_sigma8correction_positive(self, simple_pk):
        assert simple_pk.sigma8correction > 0

    def test_sigma8_positive(self, simple_pk):
        assert simple_pk.sigma8 > 0

    def test_growth_factor_normalised_at_z0(self, simple_pk):
        """D(z)/D(0) at z=0 must equal 1 by construction."""
        assert simple_pk.D(0.0) == pytest.approx(1.0, rel=1e-5)


# ─────────────────────────────────────────────────────────────────────────────
# Pz — linear power spectrum
# ─────────────────────────────────────────────────────────────────────────────

class TestPz:

    def test_positive_at_z0(self, simple_pk):
        assert np.all(simple_pk.Pz(KK, 0.0) > 0)

    def test_positive_at_z1(self, simple_pk):
        assert np.all(simple_pk.Pz(KK, 1.0) > 0)

    def test_output_shape_1d(self, simple_pk):
        result = simple_pk.Pz(KK, 0.0)
        assert result.shape == KK.shape

    def test_output_shape_2d(self, simple_pk):
        """Pz(k-array, z-array) must return shape (Nk, Nz)."""
        kk = np.array([0.01, 0.1, 1.0])
        zz = np.array([0.0, 0.5, 1.0])
        result = simple_pk.Pz(kk, zz)
        assert result.shape == (3, 3)

    def test_redshift_evolution_consistent_with_growth(self, simple_pk):
        """Pz(k,z1)/Pz(k,z0) = [D(z1)/D(z0)]^2 for all k."""
        cosmo = simple_pk.cosmo
        ratio_D2 = (cosmo.D(1.0) / cosmo.D(0.0)) ** 2
        ratio_Pk = simple_pk.Pz(KK, 1.0) / simple_pk.Pz(KK, 0.0)
        np.testing.assert_allclose(ratio_Pk, ratio_D2, rtol=1e-3)

    def test_decreasing_with_z(self, simple_pk):
        """P(k,z) < P(k,0) for z>0 (structure growth)."""
        p0 = simple_pk.Pz(KK, 0.0)
        p1 = simple_pk.Pz(KK, 1.0)
        assert np.all(p1 < p0)

    def test_regression_z0(self, simple_pk):
        expected = np.array([17.05531376, 170.55313756, 1705.53137564, 17055.31375645])
        np.testing.assert_allclose(simple_pk.Pz(KK, 0.0), expected, rtol=1e-4)

    def test_regression_z1(self, simple_pk):
        expected = np.array([6.38409121, 63.84091209, 638.40912087, 6384.09120867])
        np.testing.assert_allclose(simple_pk.Pz(KK, 1.0), expected, rtol=1e-4)

    def test_monotone_increasing_with_k(self, simple_pk):
        """P∝k is monotonically increasing (true by construction for this fixture)."""
        assert np.all(np.diff(simple_pk.Pz(KK, 0.0)) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# sigma2M and dsigma2MdM
# ─────────────────────────────────────────────────────────────────────────────

class TestSigma2M:

    def test_positive(self, simple_pk):
        assert np.all(simple_pk.sigma2M(MM, 0.0) > 0)

    def test_output_shape(self, simple_pk):
        assert simple_pk.sigma2M(MM, 0.0).shape == MM.shape

    def test_decreasing_with_M(self, simple_pk):
        """σ²(M) is a monotonically decreasing function of mass."""
        assert np.all(np.diff(simple_pk.sigma2M(MM, 0.0)) < 0)

    def test_decreasing_with_z(self, simple_pk):
        """σ²(M,z) < σ²(M,0) for z > 0 (modes grow towards z=0)."""
        s0 = simple_pk.sigma2M(MM, 0.0)
        s1 = simple_pk.sigma2M(MM, 1.0)
        assert np.all(s1 < s0)

    def test_regression_z0(self, simple_pk):
        expected = np.array([9169.25279297, 498.66831152, 86.10923449])
        np.testing.assert_allclose(simple_pk.sigma2M(MM, 0.0), expected, rtol=1e-3)

    def test_regression_z1(self, simple_pk):
        expected = np.array([3432.20576189, 186.65994828, 32.23213684])
        np.testing.assert_allclose(simple_pk.sigma2M(MM, 1.0), expected, rtol=1e-3)

    def test_dsigma2MdM_shape(self, simple_pk):
        assert simple_pk.dsigma2MdM(MM, 0.0).shape == MM.shape

    def test_dsigma2MdM_regression(self, simple_pk):
        expected = np.array([-1.10462409e-07, -6.37137776e-10, -2.12323110e-11])
        np.testing.assert_allclose(
            simple_pk.dsigma2MdM(MM, 0.0), expected, rtol=1e-3
        )


# ─────────────────────────────────────────────────────────────────────────────
# sigma2R
# ─────────────────────────────────────────────────────────────────────────────

class TestSigma2R:

    def test_positive(self, simple_pk):
        assert np.all(simple_pk.sigma2R(RR, 0.0) > 0)

    def test_output_shape(self, simple_pk):
        assert simple_pk.sigma2R(RR, 0.0).shape == RR.shape

    def test_decreasing_with_R(self, simple_pk):
        assert np.all(np.diff(simple_pk.sigma2R(RR, 0.0)) < 0)

    def test_regression(self, simple_pk):
        expected = np.array([1.89517877e+03, 4.11466363e+00, 1.94509442e-02])
        np.testing.assert_allclose(simple_pk.sigma2R(RR, 0.0), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Xi — linear correlation function
# ─────────────────────────────────────────────────────────────────────────────

class TestXi:

    def test_finite(self, simple_pk):
        rr = np.array([5.0, 10.0, 20.0])
        assert np.all(np.isfinite(simple_pk.Xi(rr, 0.0)))

    def test_output_shape(self, simple_pk):
        rr = np.array([5.0, 10.0, 20.0])
        assert simple_pk.Xi(rr, 0.0).shape == rr.shape

    def test_decreasing_with_z(self, simple_pk):
        """ξ(r,z) < ξ(r,0) for z > 0 at large r where xi > 0."""
        rr = np.array([5.0, 10.0])
        xi0 = simple_pk.Xi(rr, 0.0)
        xi1 = simple_pk.Xi(rr, 1.0)
        assert np.all(xi1 < xi0)

    def test_regression_z0(self, simple_pk):
        rr = np.array([5.0, 10.0])
        expected = np.array([1.01738916e+01, 2.41699052e+00])
        np.testing.assert_allclose(simple_pk.Xi(rr, 0.0), expected, rtol=1e-3)

    def test_regression_z1(self, simple_pk):
        rr = np.array([5.0, 10.0])
        expected = np.array([3.80825897, 0.90472026])
        np.testing.assert_allclose(simple_pk.Xi(rr, 1.0), expected, rtol=1e-3)

    def test_ratio_z0_z1_consistent_with_growth(self, simple_pk):
        """ξ(r,z1)/ξ(r,z0) ≈ [D(z1)/D(z0)]² at large scales."""
        cosmo = simple_pk.cosmo
        rr = np.array([10.0, 20.0])
        ratio_D2 = (cosmo.D(1.0) / cosmo.D(0.0)) ** 2
        xi0 = simple_pk.Xi(rr, 0.0)
        xi1 = simple_pk.Xi(rr, 1.0)
        np.testing.assert_allclose(xi1 / xi0, ratio_D2, rtol=1e-2)
