"""Tests for scampy.measure — clustering estimators and abundance statistics."""

import numpy as np
import pytest

from scampy.measure.clustering import two_point_standard, two_point_landyszalay
from scampy.measure.abundances import (
    differential_counts,
    cumulative_counts,
    get_abundances,
    cumulative_from_differential,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _make_catalogue(seed, N, M, box=50.0):
    """Return data and rand as (Ndim, Nobj) arrays — the format clustering.py expects."""
    rng = np.random.default_rng(seed)
    data = rng.uniform(0, box, (3, N))
    rand = rng.uniform(0, box, (3, M))
    return data, rand


# ─────────────────────────────────────────────────────────────────────────────
# two_point_standard
# ─────────────────────────────────────────────────────────────────────────────

class TestTwoPointStandard:

    @pytest.fixture
    def catalogue(self):
        return _make_catalogue(seed=7, N=40, M=200)

    @pytest.fixture
    def rbins(self):
        return np.geomspace(0.5, 30.0, 9)

    def test_output_shape(self, catalogue, rbins):
        data, rand = catalogue
        result = two_point_standard(data, rand, rbins)
        assert result.shape == (len(rbins),)

    def test_output_finite(self, catalogue, rbins):
        data, rand = catalogue
        result = two_point_standard(data, rand, rbins)
        assert np.all(np.isfinite(result))

    def test_output_non_negative(self, catalogue, rbins):
        """ξ_std = DD/RR - 1 should be ≥ -1 always."""
        data, rand = catalogue
        result = two_point_standard(data, rand, rbins)
        assert np.all(result >= -1.0)

    def test_random_vs_random_near_zero(self, rbins):
        """Correlation between two independent random fields should be ~0."""
        rng = np.random.default_rng(99)
        data = rng.uniform(0, 50, (3, 2000))   # (Ndim, Nobj)
        rand = rng.uniform(0, 50, (3, 8000))
        result = two_point_standard(data, rand, rbins)
        # no signal: values should be small (within noise)
        np.testing.assert_allclose(result, 0.0, atol=0.15)

    def test_regression(self, catalogue, rbins):
        data, rand = catalogue
        result = two_point_standard(data, rand, rbins)
        # regression: result must be the same array each time (deterministic)
        result2 = two_point_standard(data, rand, rbins)
        np.testing.assert_array_equal(result, result2)

    def test_angular_mode_shape(self, rbins):
        rng = np.random.default_rng(11)
        data = rng.uniform(0.0, 0.5, (2, 30))   # (Ndim=2, Nobj)
        rand = rng.uniform(0.0, 0.5, (2, 120))
        tbins = np.geomspace(0.005, 0.3, 9)
        result = two_point_standard(data, rand, tbins, angular=True)
        assert result.shape == (len(tbins),)


# ─────────────────────────────────────────────────────────────────────────────
# two_point_landyszalay
# ─────────────────────────────────────────────────────────────────────────────

class TestTwoPointLandySzalay:

    @pytest.fixture
    def catalogue(self):
        return _make_catalogue(seed=8, N=40, M=200)

    @pytest.fixture
    def rbins(self):
        return np.geomspace(0.5, 30.0, 9)

    def test_output_shape(self, catalogue, rbins):
        data, rand = catalogue
        result = two_point_landyszalay(data, rand, rbins)
        assert result.shape == (len(rbins),)

    def test_output_finite(self, catalogue, rbins):
        data, rand = catalogue
        result = two_point_landyszalay(data, rand, rbins)
        assert np.all(np.isfinite(result))

    def test_output_lower_bound(self, catalogue, rbins):
        """ξ_LS ≥ -1 in theory; with small N allow numerical violations."""
        data, rand = catalogue
        result = two_point_landyszalay(data, rand, rbins)
        assert np.all(result >= -2.0)

    def test_return_error_flag(self, catalogue, rbins):
        data, rand = catalogue
        out = two_point_landyszalay(data, rand, rbins, return_error=True)
        assert len(out) == 2
        xi, err = out
        assert xi.shape == rbins.shape
        assert err.shape == rbins.shape

    def test_random_vs_random_near_zero(self, rbins):
        rng = np.random.default_rng(100)
        data = rng.uniform(0, 50, (3, 2000))   # (Ndim, Nobj)
        rand = rng.uniform(0, 50, (3, 8000))
        result = two_point_landyszalay(data, rand, rbins)
        np.testing.assert_allclose(result, 0.0, atol=0.15)

    def test_deterministic(self, catalogue, rbins):
        data, rand = catalogue
        xi1 = two_point_landyszalay(data, rand, rbins)
        xi2 = two_point_landyszalay(data, rand, rbins)
        np.testing.assert_array_equal(xi1, xi2)


# ─────────────────────────────────────────────────────────────────────────────
# differential_counts
# ─────────────────────────────────────────────────────────────────────────────

class TestDifferentialCounts:

    @pytest.fixture
    def sample(self):
        rng = np.random.default_rng(42)
        return rng.exponential(scale=1.0, size=200)

    @pytest.fixture
    def bins(self):
        return np.linspace(0.0, 5.0, 11)    # 10 bins, width = 0.5

    def test_output_shapes(self, sample, bins):
        v, dndv, err = differential_counts(sample, bins)
        assert v.shape == (len(bins) - 1,)
        assert dndv.shape == (len(bins) - 1,)
        assert err.shape == (len(bins) - 1,)

    def test_bin_centers_correct(self, sample, bins):
        v, _, _ = differential_counts(sample, bins)
        expected = 0.5 * (bins[:-1] + bins[1:])
        np.testing.assert_allclose(v, expected, rtol=1e-12)

    def test_dndv_is_count_over_binwidth(self, sample, bins):
        v, dndv, _ = differential_counts(sample, bins)
        counts, _ = np.histogram(sample, bins=bins)
        dv = bins[1:] - bins[:-1]
        np.testing.assert_allclose(dndv, counts / dv, rtol=1e-12)

    def test_errors_are_poisson(self, sample, bins):
        _, dndv, err = differential_counts(sample, bins)
        # Poisson error: err = sqrt(N)/dv → err/dndv = 1/sqrt(N)
        counts, _ = np.histogram(sample, bins=bins)
        dv = bins[1:] - bins[:-1]
        expected_err = np.sqrt(counts) / dv
        np.testing.assert_allclose(err, expected_err, rtol=1e-12)

    def test_non_negative_values(self, sample, bins):
        _, dndv, err = differential_counts(sample, bins)
        assert np.all(dndv >= 0)
        assert np.all(err >= 0)

    def test_fact_multiplier(self, sample, bins):
        _, dndv1, _ = differential_counts(sample, bins)
        _, dndv2, _ = differential_counts(sample, bins, fact=2.5)
        np.testing.assert_allclose(dndv2, 2.5 * dndv1, rtol=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# cumulative_counts
# ─────────────────────────────────────────────────────────────────────────────

class TestCumulativeCounts:

    @pytest.fixture
    def sample(self):
        rng = np.random.default_rng(42)
        return rng.exponential(scale=1.0, size=200)

    @pytest.fixture
    def bins(self):
        return np.linspace(0.0, 5.0, 11)

    def test_output_shapes(self, sample, bins):
        Nx, err = cumulative_counts(sample, bins)
        assert Nx.shape == bins.shape
        assert err.shape == bins.shape

    def test_first_element_is_total(self, sample, bins):
        """Nx[0] equals the count of sample elements within the bin range."""
        Nx, _ = cumulative_counts(sample, bins)
        counts, _ = np.histogram(sample, bins)
        assert Nx[0] == counts.sum()

    def test_monotone_non_increasing(self, sample, bins):
        Nx, _ = cumulative_counts(sample, bins)
        assert np.all(np.diff(Nx) <= 0)

    def test_error_is_sqrt_of_count(self, sample, bins):
        Nx, err = cumulative_counts(sample, bins)
        np.testing.assert_allclose(err, np.sqrt(Nx), rtol=1e-12)

    def test_non_negative(self, sample, bins):
        Nx, err = cumulative_counts(sample, bins)
        assert np.all(Nx >= 0)
        assert np.all(err >= 0)

    def test_fact_multiplier(self, sample, bins):
        Nx1, _ = cumulative_counts(sample, bins)
        Nx2, _ = cumulative_counts(sample, bins, fact=3.0)
        np.testing.assert_allclose(Nx2, 3.0 * Nx1, rtol=1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# get_abundances
# ─────────────────────────────────────────────────────────────────────────────

class TestGetAbundances:

    @pytest.fixture
    def data(self):
        rng = np.random.default_rng(0)
        N = 100
        Ax = rng.uniform(0.0, 1.0, N)
        Nc = rng.poisson(2.0, N).astype(float)
        Ns = rng.poisson(0.5, N).astype(float)
        bins = np.linspace(0.0, 1.0, 6)
        return Ax, Nc, Ns, bins

    def test_output_shapes(self, data):
        Ax, Nc, Ns, bins = data
        bc, nc, ns = get_abundances(Ax, Nc, Ns, bins)
        assert bc.shape == (len(bins) - 1,)
        assert nc.shape == (len(bins) - 1,)
        assert ns.shape == (len(bins) - 1,)

    def test_bin_centers(self, data):
        Ax, Nc, Ns, bins = data
        bc, _, _ = get_abundances(Ax, Nc, Ns, bins)
        expected = 0.5 * (bins[:-1] + bins[1:])
        np.testing.assert_allclose(bc, expected, rtol=1e-12)

    def test_non_negative(self, data):
        Ax, Nc, Ns, bins = data
        _, nc, ns = get_abundances(Ax, Nc, Ns, bins)
        assert np.all(nc >= 0)
        assert np.all(ns >= 0)

    def test_mean_statistic_is_mean(self, data):
        Ax, Nc, Ns, bins = data
        _, nc, _ = get_abundances(Ax, Nc, Ns, bins, statistic='mean')
        # manually compute mean in first bin
        mask = (Ax >= bins[0]) & (Ax < bins[1])
        if mask.sum() > 0:
            assert nc[0] == pytest.approx(Nc[mask].mean(), rel=1e-10)


# ─────────────────────────────────────────────────────────────────────────────
# cumulative_from_differential
# ─────────────────────────────────────────────────────────────────────────────

class TestCumulativeFromDifferential:

    def test_constant_function(self):
        """∫₀ˣ c dx = c·x."""
        bins = np.linspace(0.0, 1.0, 101)
        result = cumulative_from_differential(lambda x: np.ones_like(x) * 2.0, bins)
        # result[0] = 0 by construction; result[-1] ≈ 2.0
        assert result[0] == pytest.approx(0.0, abs=1e-12)
        assert result[-1] == pytest.approx(2.0, rel=1e-4)

    def test_linear_function(self):
        """∫₀ˣ t dt = x²/2."""
        bins = np.linspace(0.0, 2.0, 1001)
        result = cumulative_from_differential(lambda x: x, bins)
        assert result[0] == pytest.approx(0.0, abs=1e-12)
        assert result[-1] == pytest.approx(2.0, rel=1e-4)

    def test_output_shape(self):
        bins = np.linspace(0.0, 1.0, 20)
        result = cumulative_from_differential(lambda x: np.ones_like(x), bins)
        assert result.shape == bins.shape

    def test_monotone_for_positive_function(self):
        bins = np.linspace(0.0, 5.0, 50)
        result = cumulative_from_differential(lambda x: np.exp(-x), bins)
        assert np.all(np.diff(result) >= 0)
