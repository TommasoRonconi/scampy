"""Tests for scampy.cosmology.model.

Covers:
- Boundary conditions (z=0 values)
- Internal consistency identities (dA = dC/(1+z), OmegaM definition, …)
- EdS analytical solutions (D ∝ a, Hz ∝ (1+z)^1.5)
- Monotonicity and sign checks
- Regression: specific values for the default ΛCDM cosmology
"""

import numpy as np
import pytest

ZZ = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 5.0])

# Reference values computed from the current implementation at ZZ above.
# Tolerances are set at rtol=1e-3 to accommodate float32 vs float64 parameter
# storage differences between the C++ and pure-Python implementations while
# still catching algorithmic errors.
REF = dict(
    Hz  = np.array([ 70.0,       73.393,   91.604,  123.248,  207.655,  566.526]),
    dC  = np.array([  0.0,      418.454, 1888.626, 3303.832, 5179.883, 7775.448]),
    dA  = np.array([  0.0,      380.413, 1259.084, 1651.916, 1726.628, 1295.908]),
    D   = np.array([  0.77898,   0.73977,   0.60229,  0.47658,  0.32828,  0.16622]),
    cdc = np.array([2.7754e11, 3.0510e11, 4.7528e11, 8.6036e11, 2.4423e12, 1.8179e13]),
    ct  = np.array([13.4503,   12.1490,    8.4097,   5.7351,    3.2100,    1.1381]),
    OM  = np.array([ 0.300,     0.363,     0.591,    0.774,     0.920,     0.989]),
    Dc  = np.array([21.99,     25.24,     36.73,    45.76,     52.83,     56.07]),
    dc  = np.array([ 1.6752,    1.6769,    1.6813,   1.6837,    1.6853,    1.6859]),
)


# ─────────────────────────────────────────────────────────────────────────────
# Boundary conditions
# ─────────────────────────────────────────────────────────────────────────────

class TestBoundaryConditions:

    def test_dC_zero_at_z0(self, default_cosmo):
        assert default_cosmo.dC(0.0) == pytest.approx(0.0, abs=1e-6)

    def test_Hz_equals_H0_at_z0(self, default_cosmo):
        assert default_cosmo.Hz(0.0) == pytest.approx(default_cosmo.H0, rel=1e-5)

    def test_H0_value(self, default_cosmo):
        assert default_cosmo.H0 == pytest.approx(70.0, rel=1e-4)

    def test_d_H0(self, default_cosmo):
        """d_H0 = c/H0 [Mpc] ≈ 4283 Mpc for h=0.7."""
        assert default_cosmo.d_H0 == pytest.approx(4282.75, rel=1e-3)

    def test_comoving_volume_zero_at_z0(self, default_cosmo):
        assert default_cosmo.comoving_volume(0.0) == pytest.approx(0.0, abs=1e-3)

    def test_comoving_volume_unit_zero_at_z0(self, default_cosmo):
        assert default_cosmo.comoving_volume_unit(0.0) == pytest.approx(0.0, abs=1e-3)

    def test_OmegaM_at_z0_equals_param(self, default_cosmo):
        assert default_cosmo.OmegaM(0.0) == pytest.approx(
            default_cosmo.param['Om_M'], rel=1e-4
        )

    def test_critical_density_comoving_at_z0(self, default_cosmo):
        """ρ_crit,comoving(0) = 2.775×10¹¹ h² M_⊙ Mpc⁻³ (within 0.1%)."""
        rho0 = default_cosmo.critical_density_comoving(0.0)
        assert rho0 == pytest.approx(2.775e11, rel=1e-3)


# ─────────────────────────────────────────────────────────────────────────────
# Internal consistency identities
# ─────────────────────────────────────────────────────────────────────────────

class TestInternalConsistency:

    def test_dA_equals_dC_over_1pz(self, default_cosmo):
        """d_A(z) = d_C(z)/(1+z) must hold exactly by construction."""
        zz = ZZ[1:]          # skip z=0 (both are 0, ratio is 0/0)
        np.testing.assert_allclose(
            default_cosmo.dA(zz),
            default_cosmo.dC(zz) / (1.0 + zz),
            rtol=1e-10,
        )

    def test_Hz_equals_H0_times_Ez(self, default_cosmo):
        """H(z) = H0 * E(z) where E²(z) = Ω_M(1+z)³ + Ω_Λ (flat)."""
        p = default_cosmo.param
        red = 1.0 + ZZ
        Ez2 = p['Om_M'] * red**3 + p['Om_L']
        expected_Hz = default_cosmo.H0 * np.sqrt(Ez2)
        np.testing.assert_allclose(default_cosmo.Hz(ZZ), expected_Hz, rtol=1e-4)

    def test_ddC_is_dH_times_zE(self, default_cosmo):
        """ddC/dz = d_H / E(z) = c / H(z), the comoving distance integrand."""
        zz = ZZ[1:]
        expected = default_cosmo.d_H0 / (default_cosmo.Hz(zz) / default_cosmo.H0)
        np.testing.assert_allclose(default_cosmo.ddC(zz), expected, rtol=1e-4)

    def test_critical_density_ratio(self, default_cosmo):
        """critical_density = critical_density_comoving * h²."""
        hh = default_cosmo.param['hh']
        np.testing.assert_allclose(
            default_cosmo.critical_density(ZZ),
            default_cosmo.critical_density_comoving(ZZ) * hh**2,
            rtol=1e-10,
        )

    def test_OmegaM_plus_OmegaL_near_one_at_high_z(self, default_cosmo):
        """At high z, matter dominates and OmegaM → 1."""
        assert default_cosmo.OmegaM(10.0) > 0.99

    def test_comoving_volume_monotone(self, default_cosmo):
        zz = np.linspace(0.01, 5.0, 50)
        vol = default_cosmo.comoving_volume(zz)
        assert np.all(np.diff(vol) > 0)

    def test_dC_monotone(self, default_cosmo):
        zz = np.linspace(0.01, 5.0, 50)
        assert np.all(np.diff(default_cosmo.dC(zz)) > 0)

    def test_Hz_monotone(self, default_cosmo):
        zz = np.linspace(0.0, 5.0, 50)
        assert np.all(np.diff(default_cosmo.Hz(zz)) > 0)

    def test_D_monotone_decreasing(self, default_cosmo):
        """Growth factor D(z) decreases with z (structure grows toward z=0)."""
        zz = np.linspace(0.0, 5.0, 50)
        assert np.all(np.diff(default_cosmo.D(zz)) < 0)

    def test_cosmic_time_monotone_decreasing(self, default_cosmo):
        """The cosmic time (age) decreases with redshift."""
        zz = np.linspace(0.01, 5.0, 50)
        assert np.all(np.diff(default_cosmo.cosmic_time(zz)) < 0)

    def test_deltac_weakly_increasing_with_z(self, default_cosmo):
        """δ_c(z) increases slightly with z (matter-dominated limit = 1.686)."""
        zz = np.array([0.0, 1.0, 5.0])
        dc = default_cosmo.deltac(zz)
        assert np.all(np.diff(dc) > 0)
        assert dc[-1] == pytest.approx(1.686, rel=5e-3)

    def test_Delta_c_monotone_increasing(self, default_cosmo):
        zz = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
        Dc = default_cosmo.Delta_c(zz)
        assert np.all(np.diff(Dc) > 0)


# ─────────────────────────────────────────────────────────────────────────────
# Einstein–de Sitter analytical solutions
# ─────────────────────────────────────────────────────────────────────────────

class TestEdSCosmology:

    def test_Hz_power_law(self, eds_cosmo):
        """EdS: H(z) = H0 (1+z)^{3/2}."""
        zz = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
        np.testing.assert_allclose(
            eds_cosmo.Hz(zz) / eds_cosmo.H0,
            (1.0 + zz) ** 1.5,
            rtol=1e-3,
        )

    def test_growth_factor_ratios(self, eds_cosmo):
        """EdS: D(z) ∝ 1/(1+z) → D(z₀)/D(z₁) = (1+z₁)/(1+z₀)."""
        zz = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
        D = eds_cosmo.D(zz)
        for i in range(1, len(zz)):
            expected_ratio = (1.0 + zz[i]) / (1.0 + zz[0])
            assert D[0] / D[i] == pytest.approx(expected_ratio, rel=5e-3)

    def test_OmegaM_equals_one_at_z0(self, eds_cosmo):
        """EdS is matter-only: Ω_M(0) = 1."""
        assert eds_cosmo.OmegaM(0.0) == pytest.approx(1.0, rel=1e-4)

    def test_OmegaM_stays_one(self, eds_cosmo):
        """EdS: Ω_M(z) = 1 for all z."""
        zz = np.array([0.0, 1.0, 5.0])
        np.testing.assert_allclose(eds_cosmo.OmegaM(zz), 1.0, rtol=1e-4)


# ─────────────────────────────────────────────────────────────────────────────
# Regression: reference values for default ΛCDM cosmology
# ─────────────────────────────────────────────────────────────────────────────

class TestRegressionDefaultCosmology:

    def test_Hz(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.Hz(ZZ), REF['Hz'], rtol=1e-3)

    def test_dC(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.dC(ZZ), REF['dC'], rtol=1e-3)

    def test_dA(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.dA(ZZ), REF['dA'], rtol=1e-3)

    def test_growth_factor(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.D(ZZ), REF['D'], rtol=1e-3)

    def test_critical_density_comoving(self, default_cosmo):
        np.testing.assert_allclose(
            default_cosmo.critical_density_comoving(ZZ), REF['cdc'], rtol=1e-3
        )

    def test_cosmic_time(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.cosmic_time(ZZ), REF['ct'], rtol=1e-3)

    def test_OmegaM(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.OmegaM(ZZ), REF['OM'], rtol=1e-2)

    def test_Delta_c(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.Delta_c(ZZ), REF['Dc'], rtol=1e-2)

    def test_deltac(self, default_cosmo):
        np.testing.assert_allclose(default_cosmo.deltac(ZZ), REF['dc'], rtol=1e-3)

    def test_param_dict_keys(self, default_cosmo):
        expected_keys = {
            'Om_M', 'Om_b', 'Om_L', 'Om_n', 'Om_r', 'Om_K',
            'hh', 'sigma8', 'w_0', 'w_a', 'zmin', 'zmax', 'thin',
        }
        assert set(default_cosmo.param.keys()) == expected_keys

    def test_zmin_zmax_attributes(self, default_cosmo):
        assert default_cosmo.zmin == pytest.approx(0.0)
        assert default_cosmo.zmax == pytest.approx(100.0, rel=1e-4)
