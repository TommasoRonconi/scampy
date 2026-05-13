"""Piecewise-linear interpolators on pre-tabulated grids.

Provides two interpolator classes sharing a common interface.

Available classes
-----------------

.. list-table::
   :widths: 30 70

   * - :class:`lin_interp`
     - Piecewise-linear interpolation in (x, y) space.
   * - :class:`log_interp`
     - Piecewise-linear interpolation in :math:`(\\log_{10} x,\\,\\log_{10} y)` space.

Both classes expose the same methods: ``__call__``, ``get_x``, ``get_y``,
and ``integrate``.
"""

import numpy


class lin_interp:
    """Piecewise-linear interpolator built from two equal-length arrays.

    Parameters
    ----------
    x : array-like
        Strictly increasing x-axis values.
    y : array-like
        Corresponding y-axis values (same length as *x*).
    """

    def __init__(self, x, y):
        self._x = numpy.asarray(x, dtype=float)
        self._y = numpy.asarray(y, dtype=float)

    def __call__(self, x):
        """Evaluate the interpolator at *x* (vectorised).

        Parameters
        ----------
        x : float or array-like
            Query point(s).

        Returns
        -------
        float or ndarray
            Interpolated value(s).  Points outside the grid are extrapolated
            linearly using the slope of the first or last segment.
        """
        scalar = numpy.ndim(x) == 0
        x = numpy.atleast_1d(numpy.asarray(x, dtype=float))
        y = numpy.interp(x, self._x, self._y)
        lo = x < self._x[0]
        if numpy.any(lo):
            slope = (self._y[1] - self._y[0]) / (self._x[1] - self._x[0])
            y[lo] = self._y[0] + slope * (x[lo] - self._x[0])
        hi = x > self._x[-1]
        if numpy.any(hi):
            slope = (self._y[-1] - self._y[-2]) / (self._x[-1] - self._x[-2])
            y[hi] = self._y[-1] + slope * (x[hi] - self._x[-1])
        return float(y[0]) if scalar else y

    def get_x(self):
        """Return the x-axis array."""
        return self._x.copy()

    def get_y(self):
        """Return the y-axis array."""
        return self._y.copy()

    def integrate(self, aa, bb):
        """Integrate the interpolated function over [aa, bb] using the trapezoidal rule.

        Parameters
        ----------
        aa : float
            Lower integration limit.
        bb : float
            Upper integration limit.

        Returns
        -------
        float
            Approximate integral.
        """
        mask = (self._x > aa) & (self._x < bb)
        xv = numpy.concatenate(([aa], self._x[mask], [bb]))
        yv = numpy.interp(xv, self._x, self._y)
        return numpy.trapz(yv, xv)


class log_interp:
    """Log-space piecewise-linear interpolator built from two equal-length arrays.

    Interpolation is performed in log10(x) vs log10(y) space.

    Parameters
    ----------
    x : array-like
        Strictly increasing positive x-axis values.
    y : array-like
        Corresponding positive y-axis values (same length as *x*).
    """

    def __init__(self, x, y):
        self._x  = numpy.asarray(x, dtype=float)
        self._y  = numpy.asarray(y, dtype=float)
        self._lx = numpy.log10(self._x)
        self._ly = numpy.log10(self._y)

    def __call__(self, x):
        """Evaluate the interpolator at *x* (vectorised).

        Parameters
        ----------
        x : float or array-like
            Query point(s).

        Returns
        -------
        float or ndarray
            Interpolated value(s).
        """
        return 10.0 ** numpy.interp(numpy.log10(x), self._lx, self._ly)

    def get_x(self):
        """Return the x-axis array."""
        return self._x.copy()

    def get_y(self):
        """Return the y-axis array."""
        return self._y.copy()

    def integrate(self, aa, bb):
        """Integrate the interpolated function over [aa, bb] using the trapezoidal rule.

        Parameters
        ----------
        aa : float
            Lower integration limit.
        bb : float
            Upper integration limit.

        Returns
        -------
        float
            Approximate integral.
        """
        mask = (self._x > aa) & (self._x < bb)
        xv = numpy.concatenate(([aa], self._x[mask], [bb]))
        yv = self(xv)
        return numpy.trapz(yv, xv)
