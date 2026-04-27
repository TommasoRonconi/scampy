"""Tests for scampy.measure.clustering_core — pair-counting kernels.

Reference values were computed with the C++ implementation using
numpy.random.default_rng(42) for the random catalogues.

Note on dA2D_DD vs dA2D_DD_omp
--------------------------------
The C++ non-OMP dA2D_DD has a known bug (the separation variable is never
assigned in the non-parallel branch).  The pure-Python implementation fixes
this; both dA2D_DD and dA2D_DD_omp are correct and equivalent there.
Tests for the angular pair counter use the reference values from dA2D_DD_omp
(which is correct in both implementations).
"""

import numpy as np
import pytest

from scampy.measure.clustering_core import (
    d2D_DD, d2D_DD_omp,
    d2D_DR, d2D_DR_omp,
    d3D_DD, d3D_DD_omp,
    d3D_DR, d3D_DR_omp,
    dA2D_DD_omp,
    dA2D_DR, dA2D_DR_omp,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers / shared data
# ─────────────────────────────────────────────────────────────────────────────

def _rng_catalogue(seed=42, N=30, M=60, box=50.0):
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, box, N)
    Y = rng.uniform(0, box, N)
    Z = rng.uniform(0, box, N)
    XR = rng.uniform(0, box, M)
    YR = rng.uniform(0, box, M)
    ZR = rng.uniform(0, box, M)
    return X, Y, Z, XR, YR, ZR


def _rng_angular(seed=42, N=30, M=60, extent=0.5):
    rng = np.random.default_rng(seed)
    RA  = rng.uniform(0, extent, N)
    Dec = rng.uniform(0, extent, N)
    RAR  = rng.uniform(0, extent, M)
    DecR = rng.uniform(0, extent, M)
    return RA, Dec, RAR, DecR


RBIN  = np.geomspace(1.0, 40.0, 9)
TBIN  = np.geomspace(0.003, 0.4, 9)

# Reference pair counts (C++ implementation, seed=42, drawn in X/Y/Z/XR/YR/ZR order)
REF_DD2  = np.array([  0,   3,   4,   6,  15,  29,  61, 110, 148])
REF_DR2  = np.array([  3,   6,  12,  23,  59, 122, 281, 408, 628])
REF_DD3  = np.array([  0,   0,   0,   1,   2,  16,  33,  99, 179])
REF_DR3  = np.array([  0,   1,   1,   5,   9,  40, 161, 349, 727])
REF_DANG = np.array([  0,   0,   0,   4,   7,  17,  49, 117, 189])
REF_DANG_DR = np.array([0,  0,   4,  10,  24,  78, 214, 494, 762])


# ─────────────────────────────────────────────────────────────────────────────
# Analytical single-pair tests
# ─────────────────────────────────────────────────────────────────────────────

class TestAnalyticalPairs:

    def test_2D_one_pair_in_correct_bin(self):
        """Two points at separation 7.5; should land in one specific bin."""
        X = np.array([0.0, 7.5])
        Y = np.array([0.0, 0.0])
        rbin = np.array([1.0, 5.0, 10.0, 20.0])
        counts = np.asarray(d2D_DD(X, Y, rbin))
        assert counts.sum() == 1
        # log10(7.5/1) / (log10(20)/3) ≈ 0.875/0.434 ≈ 2.01 → bin 2
        assert counts[2] == 1

    def test_3D_one_pair_in_correct_bin(self):
        """Two points at separation 7.5 along z-axis."""
        X = np.array([0.0, 0.0])
        Y = np.array([0.0, 0.0])
        Z = np.array([0.0, 7.5])
        rbin = np.array([1.0, 5.0, 10.0, 20.0])
        counts = np.asarray(d3D_DD(X, Y, Z, rbin))
        assert counts.sum() == 1
        assert counts[2] == 1

    def test_2D_no_pairs_outside_range(self):
        """Pairs outside [rmin, rmax] must not be counted."""
        X = np.array([0.0, 0.5])   # distance = 0.5 < rmin = 1
        Y = np.array([0.0, 0.0])
        rbin = np.array([1.0, 5.0, 10.0])
        assert np.sum(d2D_DD(X, Y, rbin)) == 0

    def test_2D_DR_one_pair(self):
        X1 = np.array([0.0])
        Y1 = np.array([0.0])
        X2 = np.array([7.5])
        Y2 = np.array([0.0])
        rbin = np.array([1.0, 5.0, 10.0, 20.0])
        counts = np.asarray(d2D_DR(X1, Y1, X2, Y2, rbin))
        assert counts.sum() == 1
        assert counts[2] == 1

    def test_3D_DR_one_pair(self):
        X1 = np.array([0.0])
        Y1 = np.array([0.0])
        Z1 = np.array([0.0])
        X2 = np.array([7.5])
        Y2 = np.array([0.0])
        Z2 = np.array([0.0])
        rbin = np.array([1.0, 5.0, 10.0, 20.0])
        counts = np.asarray(d3D_DR(X1, Y1, Z1, X2, Y2, Z2, rbin))
        assert counts.sum() == 1
        assert counts[2] == 1

    def test_total_2D_pairs_known_triangle(self):
        """Three points with separations 3, 4, 5 — all in the range [1, 10]."""
        X = np.array([0.0, 3.0, 0.0])
        Y = np.array([0.0, 0.0, 4.0])
        rbin = np.array([1.0, 10.0])
        assert np.sum(d2D_DD(X, Y, rbin)) == 3

    def test_N_choose_2_total_pairs(self):
        """For N points all within range, total pairs = N(N-1)/2."""
        rng = np.random.default_rng(0)
        N = 10
        X = rng.uniform(0, 5, N)
        Y = rng.uniform(0, 5, N)
        rbin = np.array([0.001, 100.0])
        counts = np.asarray(d2D_DD(X, Y, rbin))
        assert counts.sum() == N * (N - 1) // 2

    def test_3D_N_choose_2_total_pairs(self):
        rng = np.random.default_rng(1)
        N = 10
        X = rng.uniform(0, 5, N)
        Y = rng.uniform(0, 5, N)
        Z = rng.uniform(0, 5, N)
        rbin = np.array([0.001, 100.0])
        counts = np.asarray(d3D_DD(X, Y, Z, rbin))
        assert counts.sum() == N * (N - 1) // 2

    def test_DR_total_pairs_is_NxM(self):
        """d2D_DR with N data and M randoms gives exactly N×M pairs when all within range."""
        rng = np.random.default_rng(2)
        N, M = 5, 8
        X1 = rng.uniform(0, 3, N)
        Y1 = rng.uniform(0, 3, N)
        X2 = rng.uniform(0, 3, M)
        Y2 = rng.uniform(0, 3, M)
        rbin = np.array([0.001, 100.0])
        counts = np.asarray(d2D_DR(X1, Y1, X2, Y2, rbin))
        assert counts.sum() == N * M


# ─────────────────────────────────────────────────────────────────────────────
# OMP aliases are equivalent to serial functions
# ─────────────────────────────────────────────────────────────────────────────

class TestOmpAliases:

    def test_d2D_DD_omp_matches_serial(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(
            d2D_DD(X, Y, RBIN), d2D_DD_omp(X, Y, RBIN)
        )

    def test_d3D_DD_omp_matches_serial(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(
            d3D_DD(X, Y, Z, RBIN), d3D_DD_omp(X, Y, Z, RBIN)
        )

    def test_d3D_DR_omp_matches_serial(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(
            d3D_DR(X, Y, Z, XR, YR, ZR, RBIN),
            d3D_DR_omp(X, Y, Z, XR, YR, ZR, RBIN),
        )

    def test_dA2D_DR_omp_matches_serial(self):
        RA, Dec, RAR, DecR = _rng_angular()
        np.testing.assert_array_equal(
            dA2D_DR(RA, Dec, RAR, DecR, TBIN),
            dA2D_DR_omp(RA, Dec, RAR, DecR, TBIN),
        )


# ─────────────────────────────────────────────────────────────────────────────
# Regression against C++ reference values
# ─────────────────────────────────────────────────────────────────────────────

class TestRegressionCounts:

    def test_d2D_DD_regression(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(d2D_DD(X, Y, RBIN), REF_DD2)

    def test_d2D_DR_regression(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(d2D_DR(X, Y, XR, YR, RBIN), REF_DR2)

    def test_d3D_DD_regression(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(d3D_DD(X, Y, Z, RBIN), REF_DD3)

    def test_d3D_DR_regression(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        np.testing.assert_array_equal(d3D_DR(X, Y, Z, XR, YR, ZR, RBIN), REF_DR3)

    def test_dA2D_DD_regression(self):
        """Uses dA2D_DD_omp reference (correct in both C++ and pure-Python)."""
        RA, Dec, RAR, DecR = _rng_angular()
        np.testing.assert_array_equal(dA2D_DD_omp(RA, Dec, TBIN), REF_DANG)

    def test_dA2D_DR_regression(self):
        RA, Dec, RAR, DecR = _rng_angular()
        np.testing.assert_array_equal(dA2D_DR(RA, Dec, RAR, DecR, TBIN), REF_DANG_DR)

    def test_counts_are_non_negative(self):
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        assert np.all(np.asarray(d3D_DD(X, Y, Z, RBIN)) >= 0)
        assert np.all(np.asarray(d3D_DR(X, Y, Z, XR, YR, ZR, RBIN)) >= 0)

    def test_total_count_grows_with_range(self):
        """Pair counts for a wider bin range must be >= narrow range."""
        X, Y, Z, XR, YR, ZR = _rng_catalogue()
        rbin_narrow = np.geomspace(2.0, 20.0, 5)
        rbin_wide   = np.geomspace(1.0, 40.0, 9)
        assert (
            np.sum(d2D_DD(X, Y, rbin_wide))
            >= np.sum(d2D_DD(X, Y, rbin_narrow))
        )
