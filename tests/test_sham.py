"""Tests for scampy.sham — scatter_array, matching_sample, match_cross.

Reference values computed from the heritage C++ implementation.
No cosmology or power spectrum fixture needed.
"""

import warnings

import numpy as np
import pytest

from scampy.sham import scatter_array, matching_sample, match_cross


# ─────────────────────────────────────────────────────────────────────────────
# scatter_array
# ─────────────────────────────────────────────────────────────────────────────

class TestScatterArray:

    @pytest.fixture
    def arr(self):
        return np.array([1.0, 2.0, 3.0, 4.0, 5.0])

    def test_output_shape(self, arr):
        result = scatter_array(arr, method='normal', scale=0.5, kw_rng={'seed': 0})
        assert result.shape == arr.shape

    def test_normal_scatter_regression(self, arr):
        """Regression: fixed seed → reproducible output."""
        expected = np.array([1.15235854, 1.48000795, 3.3752256, 4.47028236, 4.02448241])
        result = scatter_array(arr, method='normal', scale=0.5, kw_rng={'seed': 42})
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_log_scatter_regression(self, arr):
        expected = np.array([1.07268377, 1.57409734, 3.56587713, 4.96725508, 3.19055684])
        result = scatter_array(arr, method='normal', log=True, scale=0.1, kw_rng={'seed': 42})
        np.testing.assert_allclose(result, expected, rtol=1e-5)

    def test_zero_scale_identity(self, arr):
        """scale=0 with normal scatter must return the original values."""
        result = scatter_array(arr, method='normal', scale=0.0, kw_rng={'seed': 1})
        np.testing.assert_allclose(result, arr, rtol=1e-12)

    def test_log_output_positive(self, arr):
        """Log-mode scatter must return strictly positive values."""
        result = scatter_array(arr, method='normal', log=True, scale=0.5, kw_rng={'seed': 5})
        assert np.all(result > 0)

    def test_normal_scatter_changes_values(self, arr):
        """With non-zero scale the output must differ from input."""
        result = scatter_array(arr, method='normal', scale=1.0, kw_rng={'seed': 99})
        assert not np.allclose(result, arr)

    def test_rng_object_accepted(self, arr):
        """Passing a Generator instance directly must work."""
        rng = np.random.default_rng(7)
        result = scatter_array(arr, method='normal', scale=0.3, rng=rng)
        assert result.shape == arr.shape

    def test_uniform_method_output_shape(self, arr):
        result = scatter_array(arr, method='uniform', scale=0.1, kw_rng={'seed': 0})
        assert result.shape == arr.shape


# ─────────────────────────────────────────────────────────────────────────────
# matching_sample
# ─────────────────────────────────────────────────────────────────────────────

class TestMatchingSample:

    def test_output_dtype_bool(self):
        mask = matching_sample(20, 5, rng=np.random.default_rng(0))
        assert mask.dtype == bool

    def test_output_size(self):
        mask = matching_sample(20, 5, rng=np.random.default_rng(0))
        assert mask.size == 20

    def test_correct_sum(self):
        """Exactly OutSize elements must be True."""
        for out in [0, 3, 7, 10]:
            mask = matching_sample(10, out, rng=np.random.default_rng(out))
            assert mask.sum() == out

    def test_regression(self):
        """Fixed seed → reproducible mask."""
        mask = matching_sample(10, 4, rng=np.random.default_rng(0))
        expected = np.array([False, False, True, False, True, True, False, True, False, False])
        np.testing.assert_array_equal(mask, expected)

    def test_full_sample(self):
        """OutSize == InSize must select all elements."""
        mask = matching_sample(8, 8, rng=np.random.default_rng(0))
        assert mask.all()

    def test_empty_sample(self):
        """OutSize == 0 must select no elements."""
        mask = matching_sample(8, 0, rng=np.random.default_rng(0))
        assert not mask.any()

    def test_weighted_sample_sums_correctly(self):
        """Even with a probability vector the sum must equal OutSize."""
        prob = np.array([0.5, 0.3, 0.1, 0.05, 0.05])
        mask = matching_sample(5, 3, InProb=prob, rng=np.random.default_rng(1))
        assert mask.sum() == 3


# ─────────────────────────────────────────────────────────────────────────────
# match_cross
# ─────────────────────────────────────────────────────────────────────────────

class TestMatchCross:

    @pytest.fixture
    def xy(self):
        rng = np.random.default_rng(7)
        X = rng.uniform(1e10, 1e14, 50)
        Y = rng.uniform(1e10, 1e14, 50)
        return X, Y

    def _run(self, X, Y, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return match_cross(X, Y, logX=True, logY=True,
                               kw_rng={'seed': 7}, **kwargs)

    def test_output_lengths_equal(self, xy):
        X, Y = xy
        ix, iy = self._run(X, Y)
        assert len(ix) == len(iy)

    def test_indices_in_range(self, xy):
        X, Y = xy
        ix, iy = self._run(X, Y)
        assert np.all((ix >= 0) & (ix < len(X)))
        assert np.all((iy >= 0) & (iy < len(Y)))

    def test_indices_unique(self, xy):
        X, Y = xy
        ix, iy = self._run(X, Y)
        assert len(np.unique(ix)) == len(ix)
        assert len(np.unique(iy)) == len(iy)

    def test_regression_first_indices(self, xy):
        X, Y = xy
        ix, iy = self._run(X, Y)
        np.testing.assert_array_equal(ix[:5], [16, 19, 27,  1, 41])
        np.testing.assert_array_equal(iy[:5], [ 8, 18, 26, 31, 27])

    def test_regression_lengths(self, xy):
        X, Y = xy
        ix, iy = self._run(X, Y)
        assert len(ix) == 50
        assert len(iy) == 50

    def test_deterministic(self, xy):
        """Two calls with the same seed must return identical results."""
        X, Y = xy
        ix1, iy1 = self._run(X, Y)
        ix2, iy2 = self._run(X, Y)
        np.testing.assert_array_equal(ix1, ix2)
        np.testing.assert_array_equal(iy1, iy2)

    def test_linear_scale(self, xy):
        """match_cross with logX=False, logY=False must not raise."""
        X, Y = xy
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ix, iy = match_cross(X, Y, kw_rng={'seed': 0})
        assert len(ix) == len(iy)
