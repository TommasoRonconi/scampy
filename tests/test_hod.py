"""Tests for scampy.hod — HOD, HOD_unconditioned_sat, HOD_zdep.

Reference values computed from the heritage C++ implementation.
No power spectrum or cosmology fixture is required; HOD is pure math.
"""

import numpy as np
import pytest

from scampy.hod import HOD, HOD_unconditioned_sat, HOD_zdep


MH = np.array([1e11, 1e12, 1e13, 1e14])


# ─────────────────────────────────────────────────────────────────────────────
# HOD
# ─────────────────────────────────────────────────────────────────────────────

class TestHODPcen:

    @pytest.fixture
    def hod(self):
        return HOD(Mmin=1e12, sigma=0.5, M0=5e11, M1=1e13, alpha=1.0)

    def test_output_shape(self, hod):
        assert hod.Pcen(MH).shape == MH.shape

    def test_range_zero_to_one(self, hod):
        assert np.all(hod.Pcen(MH) >= 0.0)
        assert np.all(hod.Pcen(MH) <= 1.0)

    def test_monotone_increasing(self, hod):
        """Pcen is a sigmoid — monotonically increasing with Mh."""
        assert np.all(np.diff(hod.Pcen(MH)) > 0)

    def test_at_mmin_equals_half(self, hod):
        """By definition Pcen(Mmin) = 0.5."""
        result = hod.Pcen(np.array([1e12]))
        assert result.item() == pytest.approx(0.5, rel=1e-6)

    def test_below_mmin_near_zero(self, hod):
        """Pcen ≪ 1 for Mh ≪ Mmin."""
        result = hod.Pcen(np.array([1e9, 1e10, 1e11]))
        assert np.all(result < 0.1)

    def test_above_mmin_near_one(self, hod):
        """Pcen → 1 for Mh ≫ Mmin."""
        result = hod.Pcen(np.array([1e14, 1e15]))
        assert np.all(result > 0.99)

    def test_zero_mass_returns_zero(self, hod):
        assert hod.Pcen(np.array([0.0]))[0] == 0.0

    def test_regression(self, hod):
        expected = np.array([0.00233887, 0.5, 0.99766113, 0.99999999])
        np.testing.assert_allclose(hod.Pcen(MH), expected, rtol=1e-4)


class TestHODPsat:

    @pytest.fixture
    def hod(self):
        return HOD(Mmin=1e12, sigma=0.5, M0=5e11, M1=1e13, alpha=1.0)

    def test_output_shape(self, hod):
        assert hod.Psat(MH).shape == MH.shape

    def test_non_negative(self, hod):
        assert np.all(hod.Psat(MH) >= 0.0)

    def test_zero_below_M0(self, hod):
        """Psat = 0 for all Mh ≤ M0."""
        hod2 = HOD(Mmin=1e12, M0=1e13, M1=1e14, alpha=1.0)
        result = hod2.Psat(np.array([1e11, 1e12, 1e13]))
        np.testing.assert_array_equal(result, 0.0)

    def test_zero_mass_returns_zero(self, hod):
        assert hod.Psat(np.array([0.0]))[0] == 0.0

    def test_conditioned_on_pcen(self, hod):
        """Psat ≤ Pcen × (something) — conditioned HOD has Psat ≤ Pcen × max_factor."""
        # At Mh ≫ M0 and Mh ≫ Mmin, Psat = Pcen × ((Mh-M0)/M1)^alpha
        # At Mh = 1e14: Pcen≈1, Psat ≈ ((1e14-5e11)/1e13)^1 ≈ 9.95
        assert hod.Psat(np.array([1e14]))[0] == pytest.approx(9.94999992, rel=1e-4)

    def test_regression(self, hod):
        expected = np.array([0.0, 0.025, 0.94777808, 9.94999992])
        np.testing.assert_allclose(hod.Psat(MH), expected, rtol=1e-4)

    def test_increasing_with_mass_above_M0(self, hod):
        """Psat increases with Mh where Mh > M0."""
        masses = np.geomspace(1e12, 1e15, 10)
        psat = hod.Psat(masses)
        assert np.all(np.diff(psat) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# HOD_unconditioned_sat
# ─────────────────────────────────────────────────────────────────────────────

class TestHODUnconditionedSat:

    @pytest.fixture
    def hod(self):
        return HOD_unconditioned_sat(Mmin=1e12, sigma=0.5, M0=5e11, M1=1e13, alpha=1.0)

    @pytest.fixture
    def hod_cond(self):
        return HOD(Mmin=1e12, sigma=0.5, M0=5e11, M1=1e13, alpha=1.0)

    def test_pcen_identical_to_hod(self, hod, hod_cond):
        """Unconditioned and conditioned HOD share the same Pcen."""
        np.testing.assert_allclose(hod.Pcen(MH), hod_cond.Pcen(MH), rtol=1e-12)

    def test_psat_differs_from_hod(self, hod, hod_cond):
        """Unconditioned Psat is NOT multiplied by Pcen — differs from HOD."""
        # Unconditioned: (Mh-M0)/M1; Conditioned: Pcen * (Mh-M0)/M1
        # At Mh=1e11 where Pcen < 0.5, the conditioned value is smaller.
        # At Mh=1e14 where Pcen≈1, values are nearly equal.
        # Below Mmin, unconditioned > conditioned.
        psat_uc  = hod.Psat(np.array([1e13]))
        psat_con = hod_cond.Psat(np.array([1e13]))
        assert psat_uc[0] > psat_con[0]

    def test_psat_non_negative(self, hod):
        assert np.all(hod.Psat(MH) >= 0.0)

    def test_psat_regression(self, hod):
        expected = np.array([0.0, 0.05, 0.95, 9.95])
        np.testing.assert_allclose(hod.Psat(MH), expected, rtol=1e-4)


# ─────────────────────────────────────────────────────────────────────────────
# HOD_zdep
# ─────────────────────────────────────────────────────────────────────────────

class TestHODZdep:

    @pytest.fixture
    def hz(self):
        return HOD_zdep(
            redshift=0.5,
            lMmin=12.0, lsigma=-0.3,
            lM0=11.7, lM1=13.0, lalpha=0.0
        )

    def test_pcen_shape(self, hz):
        assert hz.Pcen(MH).shape == MH.shape

    def test_pcen_range(self, hz):
        result = hz.Pcen(MH)
        assert np.all(result >= 0.0)
        assert np.all(result <= 1.0)

    def test_pcen_at_mmin(self, hz):
        """Pcen(10^lMmin) = 0.5 regardless of redshift."""
        Mmin = 10 ** 12.0
        result = hz.Pcen(np.array([Mmin]))
        assert result.item() == pytest.approx(0.5, rel=1e-4)

    def test_psat_non_negative(self, hz):
        assert np.all(hz.Psat(MH) >= 0.0)

    def test_pcen_regression_z05(self, hz):
        expected = np.array([0.00238829, 0.5, 0.99761171, 0.99999999])
        np.testing.assert_allclose(hz.Pcen(MH), expected, rtol=1e-4)

    def test_psat_regression_z05(self, hz):
        expected = np.array([0.0, 0.02494064, 0.94761268, 9.94988119])
        np.testing.assert_allclose(hz.Psat(MH), expected, rtol=1e-4)

    def test_set_changes_redshift(self, hz):
        """Calling set(redshift=z) updates the internal redshift."""
        hz.set(redshift=1.0)
        assert hz.z[0] == pytest.approx(1.0, rel=1e-12)

    def test_zero_slopes_redshift_independent(self):
        """With all _b=0 slopes, Pcen/Psat must not change with redshift."""
        hz0 = HOD_zdep(redshift=0.0, lMmin=12.0, lsigma=0.0)
        hz1 = HOD_zdep(redshift=1.0, lMmin=12.0, lsigma=0.0)
        np.testing.assert_allclose(hz0.Pcen(MH), hz1.Pcen(MH), rtol=1e-8)
        np.testing.assert_allclose(hz0.Psat(MH), hz1.Psat(MH), rtol=1e-8)

    def test_matches_hod_at_z0_with_zero_slopes(self):
        """HOD_zdep with slopes=0 and z=0 must equal plain HOD at z=0."""
        hz = HOD_zdep(
            redshift=0.0,
            lMmin=12.0, lsigma=-0.3,
            lM0=11.7, lM1=13.0, lalpha=0.0
        )
        hod = HOD(
            Mmin=1e12, sigma=10**(-0.3),
            M0=10**11.7, M1=10**13.0, alpha=1.0
        )
        np.testing.assert_allclose(hz.Pcen(MH), hod.Pcen(MH), rtol=1e-5)
