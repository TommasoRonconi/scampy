"""Tests for scampy.halo.model.halo_model.

All tests use the simple_pk (P∝k) fixture and a standard HOD.
Reference values computed from the heritage C++ implementation with thin=128.

Note: with P∝k the 1-halo term is negligibly small at the scales tested
and the Xi/Wr functions can be negative at small separations (Fourier
artefact of the unrealistic spectrum).  Regression tests capture the
actual reference values; qualitative tests are limited to well-behaved
observables (ng, bias, Nz normalisation, Pk decomposition).
"""

import numpy as np
import pytest

from scampy.halo.model import halo_model
from scampy.hod import HOD
from scampy.power_spectrum import power_spectrum


@pytest.fixture(scope="module")
def hod():
    return HOD(Mmin=1e12, sigma=0.5, M0=5e11, M1=1e13, alpha=1.0)


@pytest.fixture(scope="module")
def hm(simple_pk):
    return halo_model(simple_pk, thin=128)


KK = np.geomspace(0.01, 10.0, 6)
RR = np.array([1.0, 5.0, 20.0])
TH = np.array([1e-3, 5e-3, 2e-2])
ZBINS = np.linspace(0.1, 1.0, 6)


# ─────────────────────────────────────────────────────────────────────────────
# Construction
# ─────────────────────────────────────────────────────────────────────────────

class TestConstruction:

    def test_wrong_pk_raises(self):
        with pytest.raises((ValueError, TypeError)):
            halo_model("not a power_spectrum")

    def test_mass_grid_shape(self, hm):
        assert hm.Mh_grid.shape == (128,)

    def test_k_grid_shape(self, hm):
        assert hm.kh_grid.shape == (128,)

    def test_mass_grid_positive(self, hm):
        assert np.all(hm.Mh_grid > 0)


# ─────────────────────────────────────────────────────────────────────────────
# ng — galaxy number density
# ─────────────────────────────────────────────────────────────────────────────

class TestNg:

    def test_positive_at_z0(self, hm, hod):
        assert hm.ng(hod, 0.0) > 0

    def test_finite(self, hm, hod):
        assert np.isfinite(hm.ng(hod, 0.0))

    def test_array_shape(self, hm, hod):
        zz = np.array([0.0, 0.5, 1.0])
        result = hm.ng(hod, zz)
        assert result.shape == zz.shape

    def test_positive_array(self, hm, hod):
        zz = np.array([0.0, 0.5, 1.0])
        assert np.all(hm.ng(hod, zz) > 0)

    def test_regression_z0(self, hm, hod):
        assert hm.ng(hod, 0.0) == pytest.approx(0.02102019861448587, rel=1e-3)

    def test_regression_z1(self, hm, hod):
        assert hm.ng(hod, 1.0) == pytest.approx(0.01950635750482269, rel=1e-3)

    def test_regression_array(self, hm, hod):
        zz = np.array([0.0, 0.5, 1.0])
        expected = np.array([0.0210202, 0.01991621, 0.01950636])
        np.testing.assert_allclose(hm.ng(hod, zz), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# bias — effective galaxy bias
# ─────────────────────────────────────────────────────────────────────────────

class TestBias:

    def test_positive_at_z0(self, hm, hod):
        assert hm.bias(hod, 0.0) > 0

    def test_finite(self, hm, hod):
        assert np.isfinite(hm.bias(hod, 0.0))

    def test_regression_z0(self, hm, hod):
        assert hm.bias(hod, 0.0) == pytest.approx(0.7894284256393715, rel=1e-3)

    def test_regression_z1(self, hm, hod):
        assert hm.bias(hod, 1.0) == pytest.approx(0.7247928875492371, rel=1e-3)

    def test_array_shape(self, hm, hod):
        zz = np.array([0.0, 0.5, 1.0])
        result = hm.bias(hod, zz)
        assert result.shape == zz.shape


# ─────────────────────────────────────────────────────────────────────────────
# Mhalo — mean host halo mass
# ─────────────────────────────────────────────────────────────────────────────

class TestMhalo:

    def test_positive(self, hm, hod):
        assert hm.Mhalo(hod, 0.0) > 0

    def test_finite(self, hm, hod):
        assert np.isfinite(hm.Mhalo(hod, 0.0))

    def test_above_mmin(self, hm, hod):
        """Mean host halo mass must exceed the HOD minimum mass."""
        assert hm.Mhalo(hod, 0.0) > hod.Mmin

    def test_regression(self, hm, hod):
        assert hm.Mhalo(hod, 0.0) == pytest.approx(20327835045261.145, rel=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# dngdM — differential galaxy number density
# ─────────────────────────────────────────────────────────────────────────────

class TestDngdM:

    def test_output_shape(self, hm, hod):
        mm = np.array([1e12, 1e13, 1e14])
        result = hm.dngdM(mm, hod, 0.0)
        assert result.shape == mm.shape

    def test_non_negative(self, hm, hod):
        mm = np.array([1e12, 1e13, 1e14])
        assert np.all(hm.dngdM(mm, hod, 0.0) >= 0)

    def test_regression(self, hm, hod):
        mm = np.array([1e12, 1e13, 1e14])
        expected = np.array([5.39751005e-15, 3.10213264e-16, 2.09597625e-17])
        np.testing.assert_allclose(hm.dngdM(mm, hod, 0.0), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Pk — galaxy power spectrum
# ─────────────────────────────────────────────────────────────────────────────

class TestPk:

    def test_Pk_equals_sum_of_terms(self, hm, hod):
        """Pk = Pk_1halo + Pk_2halo must hold exactly."""
        pk1 = hm.Pk_1halo(KK, hod, 0.0)
        pk2 = hm.Pk_2halo(KK, hod, 0.0)
        pk  = hm.Pk(KK, hod, 0.0)
        np.testing.assert_allclose(pk, pk1 + pk2, rtol=1e-10)

    def test_2halo_dominant_at_low_k(self, hm, hod):
        """At large scales (low k) the 2-halo term dominates."""
        kk_low = np.array([0.01, 0.05])
        pk1 = hm.Pk_1halo(kk_low, hod, 0.0)
        pk2 = hm.Pk_2halo(kk_low, hod, 0.0)
        assert np.all(pk2 > pk1)

    def test_Pk_2halo_positive(self, hm, hod):
        assert np.all(hm.Pk_2halo(KK, hod, 0.0) > 0)

    def test_Pk_1halo_non_negative(self, hm, hod):
        assert np.all(hm.Pk_1halo(KK, hod, 0.0) >= 0)

    def test_output_shape(self, hm, hod):
        assert hm.Pk(KK, hod, 0.0).shape == KK.shape
        assert hm.Pk_1halo(KK, hod, 0.0).shape == KK.shape
        assert hm.Pk_2halo(KK, hod, 0.0).shape == KK.shape

    def test_regression_Pk_1halo(self, hm, hod):
        expected = np.array([0., 0., 0., 0., 0., 2.93530291e-08])
        np.testing.assert_allclose(hm.Pk_1halo(KK, hod, 0.0), expected, rtol=1e-3,
                                   atol=1e-12)

    def test_regression_Pk_2halo(self, hm, hod):
        expected = np.array([10.62868726, 42.30545803, 167.91135963,
                              638.10700046, 1618.50627409, 1799.95263547])
        np.testing.assert_allclose(hm.Pk_2halo(KK, hod, 0.0), expected, rtol=1e-3)

    def test_regression_Pk(self, hm, hod):
        expected = np.array([10.62868726, 42.30545803, 167.91135963,
                              638.10700046, 1618.50627409, 1799.9526355])
        np.testing.assert_allclose(hm.Pk(KK, hod, 0.0), expected, rtol=1e-3)

    def test_Pk_at_z05_regression(self, hm, hod):
        kk = np.geomspace(0.01, 10.0, 4)
        expected = np.array([5.20289403, 52.00997562, 501.97415964, 1532.61877302])
        np.testing.assert_allclose(hm.Pk(kk, hod, 0.5), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Xi — 3D correlation function
# ─────────────────────────────────────────────────────────────────────────────

class TestXi:

    def test_output_shape(self, hm, hod):
        assert hm.Xi(RR, hod, 0.0).shape == RR.shape
        assert hm.Xi_1halo(RR, hod, 0.0).shape == RR.shape
        assert hm.Xi_2halo(RR, hod, 0.0).shape == RR.shape

    def test_finite(self, hm, hod):
        assert np.all(np.isfinite(hm.Xi(RR, hod, 0.0)))
        assert np.all(np.isfinite(hm.Xi_2halo(RR, hod, 0.0)))

    def test_Xi_equals_sum_of_terms(self, hm, hod):
        """Xi ≈ Xi_1halo + Xi_2halo (via shared FFT — may differ by ~1e-4)."""
        xi1 = hm.Xi_1halo(RR, hod, 0.0)
        xi2 = hm.Xi_2halo(RR, hod, 0.0)
        xi  = hm.Xi(RR, hod, 0.0)
        np.testing.assert_allclose(xi, xi1 + xi2, rtol=1e-3)

    def test_regression_Xi_1halo(self, hm, hod):
        expected = np.array([-6.15389482e-06, 4.83540038e-07, 4.74203802e-07])
        np.testing.assert_allclose(hm.Xi_1halo(RR, hod, 0.0), expected, rtol=1e-3,
                                   atol=1e-10)

    def test_regression_Xi_2halo(self, hm, hod):
        expected = np.array([-6.71070069e+01, -1.50816549e-01, 3.85889616e-02])
        np.testing.assert_allclose(hm.Xi_2halo(RR, hod, 0.0), expected, rtol=1e-3)

    def test_regression_Xi(self, hm, hod):
        expected = np.array([-6.71070131e+01, -1.50816066e-01, 3.85894358e-02])
        np.testing.assert_allclose(hm.Xi(RR, hod, 0.0), expected, rtol=1e-3)

    def test_z_broadcasting_raises(self, hm, hod):
        """Xi does not support z-array input; must raise RuntimeError."""
        with pytest.raises(RuntimeError):
            hm.Xi(RR, hod, np.array([0.0, 0.5]))


# ─────────────────────────────────────────────────────────────────────────────
# Wr — projected correlation function
# ─────────────────────────────────────────────────────────────────────────────

class TestWr:

    def test_output_shape(self, hm, hod):
        assert hm.Wr(RR, hod, 0.0).shape == RR.shape

    def test_finite(self, hm, hod):
        assert np.all(np.isfinite(hm.Wr(RR, hod, 0.0)))

    def test_regression_Wr_1halo(self, hm, hod):
        expected = np.array([1.16735923e-05, 9.15217558e-06, 9.15362888e-06])
        np.testing.assert_allclose(hm.Wr_1halo(RR, hod, 0.0), expected, rtol=1e-3)

    def test_regression_Wr_2halo(self, hm, hod):
        expected = np.array([-119.34275546,   -0.8252957,    0.60964273])
        np.testing.assert_allclose(hm.Wr_2halo(RR, hod, 0.0), expected, rtol=1e-3)

    def test_regression_Wr(self, hm, hod):
        expected = np.array([-119.34274379,   -0.82528655,    0.60965188])
        np.testing.assert_allclose(hm.Wr(RR, hod, 0.0), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Wt — angular correlation function
# ─────────────────────────────────────────────────────────────────────────────

class TestWt:

    def test_output_shape(self, hm, hod):
        assert hm.Wt(TH, hod, 1.0).shape == TH.shape

    def test_finite(self, hm, hod):
        assert np.all(np.isfinite(hm.Wt(TH, hod, 1.0)))

    def test_regression_Wt_1halo(self, hm, hod):
        expected = np.array([1.84658014e-05, 1.84201950e-05, 1.84231323e-05])
        np.testing.assert_allclose(hm.Wt_1halo(TH, hod, 1.0), expected, rtol=1e-3)

    def test_regression_Wt_2halo(self, hm, hod):
        expected = np.array([-1.09729902,  0.41233833,  0.42426542])
        np.testing.assert_allclose(hm.Wt_2halo(TH, hod, 1.0), expected, rtol=1e-3)

    def test_regression_Wt(self, hm, hod):
        expected = np.array([-1.09728055,  0.41235675,  0.42428384])
        np.testing.assert_allclose(hm.Wt(TH, hod, 1.0), expected, rtol=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Nz — redshift distribution
# ─────────────────────────────────────────────────────────────────────────────

class TestNz:

    def test_output_shape(self, hm, hod):
        result = hm.Nz(hod, ZBINS)
        assert result.shape == (len(ZBINS) - 1,)

    def test_non_negative(self, hm, hod):
        result = hm.Nz(hod, ZBINS)
        assert np.all(result >= 0)

    def test_normalised(self, hm, hod):
        """N(z) dz must sum to 1 (normalised redshift distribution)."""
        Nz = hm.Nz(hod, ZBINS)
        dz = np.diff(ZBINS)
        assert (Nz * dz).sum() == pytest.approx(1.0, rel=1e-6)

    def test_regression(self, hm, hod):
        expected = np.array([0.20380512, 0.62697017, 1.12000873, 1.59534808, 2.00942346])
        np.testing.assert_allclose(hm.Nz(hod, ZBINS), expected, rtol=1e-3)
