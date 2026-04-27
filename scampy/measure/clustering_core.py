"""Pair-counting kernels for two-point statistics.

Provides brute-force :math:`O(N^2)` pair counters used by the two-point estimators
in :mod:`scampy.measure.clustering`.  Pairs are binned in log-spaced
separation bins; the ``_omp`` variants are provided for API compatibility
and are identical to their serial counterparts.

Available functions
-------------------

**2D Cartesian pair counts**

.. list-table::
   :widths: 35 65

   * - :func:`d2D_DD`, :func:`d2D_DD_omp`
     - Data–data pair counts on a 2D Cartesian catalogue.
   * - :func:`d2D_DR`, :func:`d2D_DR_omp`
     - Data–random pair counts on two 2D Cartesian catalogues.

**2D angular pair counts**

.. list-table::
   :widths: 35 65

   * - :func:`dA2D_DD`, :func:`dA2D_DD_omp`
     - Data–data pair counts using angular separations.
   * - :func:`dA2D_DR`, :func:`dA2D_DR_omp`
     - Data–random pair counts using angular separations.

**3D Cartesian pair counts**

.. list-table::
   :widths: 35 65

   * - :func:`d3D_DD`, :func:`d3D_DD_omp`
     - Data–data pair counts on a 3D Cartesian catalogue.
   * - :func:`d3D_DR`, :func:`d3D_DR_omp`
     - Data–random pair counts on two 3D Cartesian catalogues.
"""

import numpy


# ---------------------------------------------------------------------------
# Internal helper

def _log_bin(rr, rmin, delta, nbins):
    """Return bin indices for separations within [rmin, rmax]."""
    ib = (numpy.log10(rr / rmin) / delta).astype(int)
    return numpy.clip(ib, 0, nbins - 1)


# ---------------------------------------------------------------------------
# 2D Cartesian

def d2D_DD(X, Y, rbin):
    """Count data-data pairs in 2D separation bins.

    Parameters
    ----------
    X : array-like of float
        X-coordinates of the catalogue.
    Y : array-like of float
        Y-coordinates of the catalogue.
    rbin : array-like of float
        Separation bin edges (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    X, Y = numpy.asarray(X, dtype=float), numpy.asarray(Y, dtype=float)
    rbin = numpy.asarray(rbin, dtype=float)
    rmin, rmax = rbin[0], rbin[-1]
    nbins = len(rbin)
    delta = numpy.log10(rmax / rmin) / nbins
    NDD = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(X)):
        dx = X[ii] - X[ii + 1:]
        dy = Y[ii] - Y[ii + 1:]
        rr = numpy.sqrt(dx * dx + dy * dy)
        mask = (rmin <= rr) & (rr <= rmax)
        if mask.any():
            NDD += numpy.bincount(_log_bin(rr[mask], rmin, delta, nbins),
                                  minlength=nbins)
    return NDD


def d2D_DR(X1, Y1, X2, Y2, rbin):
    """Count data-random cross-pairs in 2D separation bins.

    Parameters
    ----------
    X1 : array-like of float
        X-coordinates of the first (data) catalogue.
    Y1 : array-like of float
        Y-coordinates of the first catalogue.
    X2 : array-like of float
        X-coordinates of the second (random) catalogue.
    Y2 : array-like of float
        Y-coordinates of the second catalogue.
    rbin : array-like of float
        Separation bin edges (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    X1, Y1 = numpy.asarray(X1, dtype=float), numpy.asarray(Y1, dtype=float)
    X2, Y2 = numpy.asarray(X2, dtype=float), numpy.asarray(Y2, dtype=float)
    rbin = numpy.asarray(rbin, dtype=float)
    rmin, rmax = rbin[0], rbin[-1]
    nbins = len(rbin)
    delta = numpy.log10(rmax / rmin) / nbins
    NDR = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(X1)):
        dx = X1[ii] - X2
        dy = Y1[ii] - Y2
        rr = numpy.sqrt(dx * dx + dy * dy)
        mask = (rmin <= rr) & (rr <= rmax)
        if mask.any():
            NDR += numpy.bincount(_log_bin(rr[mask], rmin, delta, nbins),
                                  minlength=nbins)
    return NDR


# ---------------------------------------------------------------------------
# 2D angular

def dA2D_DD(RA, Dec, thetabin):
    """Count data-data pairs in 2D angular separation bins.

    Parameters
    ----------
    RA : array-like of float
        Right-ascension coordinates [rad].
    Dec : array-like of float
        Declination coordinates [rad].
    thetabin : array-like of float
        Angular separation bin edges [rad] (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    RA  = numpy.asarray(RA,  dtype=float)
    Dec = numpy.asarray(Dec, dtype=float)
    thetabin = numpy.asarray(thetabin, dtype=float)
    tmin, tmax = thetabin[0], thetabin[-1]
    nbins = len(thetabin)
    delta = numpy.log10(tmax / tmin) / nbins
    NDD = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(RA)):
        dra  = (RA[ii] - RA[ii + 1:]) * numpy.cos(0.5 * (Dec[ii] + Dec[ii + 1:]))
        ddec = Dec[ii] - Dec[ii + 1:]
        tt   = numpy.sqrt(dra * dra + ddec * ddec)
        mask = (tmin <= tt) & (tt <= tmax)
        if mask.any():
            NDD += numpy.bincount(_log_bin(tt[mask], tmin, delta, nbins),
                                  minlength=nbins)
    return NDD


def dA2D_DR(RA1, Dec1, RA2, Dec2, thetabin):
    """Count data-random cross-pairs in 2D angular separation bins.

    Parameters
    ----------
    RA1 : array-like of float
        Right-ascension of the first (data) catalogue [rad].
    Dec1 : array-like of float
        Declination of the first catalogue [rad].
    RA2 : array-like of float
        Right-ascension of the second (random) catalogue [rad].
    Dec2 : array-like of float
        Declination of the second catalogue [rad].
    thetabin : array-like of float
        Angular separation bin edges [rad] (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    RA1,  Dec1 = numpy.asarray(RA1,  dtype=float), numpy.asarray(Dec1, dtype=float)
    RA2,  Dec2 = numpy.asarray(RA2,  dtype=float), numpy.asarray(Dec2, dtype=float)
    thetabin = numpy.asarray(thetabin, dtype=float)
    tmin, tmax = thetabin[0], thetabin[-1]
    nbins = len(thetabin)
    delta = numpy.log10(tmax / tmin) / nbins
    NDR = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(RA1)):
        dra  = (RA1[ii] - RA2) * numpy.cos(0.5 * (Dec1[ii] + Dec2))
        ddec = Dec1[ii] - Dec2
        tt   = numpy.sqrt(dra * dra + ddec * ddec)
        mask = (tmin <= tt) & (tt <= tmax)
        if mask.any():
            NDR += numpy.bincount(_log_bin(tt[mask], tmin, delta, nbins),
                                  minlength=nbins)
    return NDR


# ---------------------------------------------------------------------------
# 3D Cartesian

def d3D_DD(X, Y, Z, rbin):
    """Count data-data pairs in 3D separation bins.

    Parameters
    ----------
    X : array-like of float
        X-coordinates of the catalogue.
    Y : array-like of float
        Y-coordinates of the catalogue.
    Z : array-like of float
        Z-coordinates of the catalogue.
    rbin : array-like of float
        Separation bin edges (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    X = numpy.asarray(X, dtype=float)
    Y = numpy.asarray(Y, dtype=float)
    Z = numpy.asarray(Z, dtype=float)
    rbin = numpy.asarray(rbin, dtype=float)
    rmin, rmax = rbin[0], rbin[-1]
    nbins = len(rbin)
    delta = numpy.log10(rmax / rmin) / nbins
    NDD = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(X)):
        dx = X[ii] - X[ii + 1:]
        dy = Y[ii] - Y[ii + 1:]
        dz = Z[ii] - Z[ii + 1:]
        rr = numpy.sqrt(dx * dx + dy * dy + dz * dz)
        mask = (rmin <= rr) & (rr <= rmax)
        if mask.any():
            NDD += numpy.bincount(_log_bin(rr[mask], rmin, delta, nbins),
                                  minlength=nbins)
    return NDD


def d3D_DR(X1, Y1, Z1, X2, Y2, Z2, rbin):
    """Count data-random cross-pairs in 3D separation bins.

    Parameters
    ----------
    X1 : array-like of float
        X-coordinates of the first (data) catalogue.
    Y1 : array-like of float
        Y-coordinates of the first catalogue.
    Z1 : array-like of float
        Z-coordinates of the first catalogue.
    X2 : array-like of float
        X-coordinates of the second (random) catalogue.
    Y2 : array-like of float
        Y-coordinates of the second catalogue.
    Z2 : array-like of float
        Z-coordinates of the second catalogue.
    rbin : array-like of float
        Separation bin edges (log-spaced).

    Returns
    -------
    ndarray of int
        Pair counts per bin.
    """
    X1 = numpy.asarray(X1, dtype=float)
    Y1 = numpy.asarray(Y1, dtype=float)
    Z1 = numpy.asarray(Z1, dtype=float)
    X2 = numpy.asarray(X2, dtype=float)
    Y2 = numpy.asarray(Y2, dtype=float)
    Z2 = numpy.asarray(Z2, dtype=float)
    rbin = numpy.asarray(rbin, dtype=float)
    rmin, rmax = rbin[0], rbin[-1]
    nbins = len(rbin)
    delta = numpy.log10(rmax / rmin) / nbins
    NDR = numpy.zeros(nbins, dtype=numpy.intp)
    for ii in range(len(X1)):
        dx = X1[ii] - X2
        dy = Y1[ii] - Y2
        dz = Z1[ii] - Z2
        rr = numpy.sqrt(dx * dx + dy * dy + dz * dz)
        mask = (rmin <= rr) & (rr <= rmax)
        if mask.any():
            NDR += numpy.bincount(_log_bin(rr[mask], rmin, delta, nbins),
                                  minlength=nbins)
    return NDR


# ---------------------------------------------------------------------------
# _omp aliases — API-compatible, serial implementation

d2D_DD_omp  = d2D_DD
d2D_DR_omp  = d2D_DR
dA2D_DD_omp = dA2D_DD
dA2D_DR_omp = dA2D_DR
d3D_DD_omp  = d3D_DD
d3D_DR_omp  = d3D_DR
