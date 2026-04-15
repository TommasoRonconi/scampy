#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <interpolation.h>
#include <vector>

namespace py = pybind11;

template class utl::interpolator< utl::lin_interp >;
template class utl::interpolator< utl::log_interp >;

#define INTERP_INIT_DOC \
  "\nParameters\n----------\n" \
  "x : list of float\n    Strictly increasing x-axis values.\n" \
  "y : list of float\n    Corresponding y-axis values (same length as x)."

#define INTEGRATE_DOC \
  "Integrate the interpolated function over [aa, bb] using the trapezoidal rule.\n" \
  "\nParameters\n----------\n" \
  "aa : float\n    Lower integration limit.\n" \
  "bb : float\n    Upper integration limit.\n" \
  "\nReturns\n-------\nfloat\n    Approximate integral."

PYBIND11_MODULE( interpolation, m ) {

  py::class_< utl::interpolator< utl::lin_interp > >( m, "lin_interp",
    "Piecewise-linear interpolator built from two equal-length arrays.\n"
    INTERP_INIT_DOC )
    .def(py::init< const std::vector< double > &, const std::vector< double > & >())
    .def("get_x", &utl::interpolator< utl::lin_interp >::get_xv,
	 "Return the x-axis array." )
    .def("get_y", &utl::interpolator< utl::lin_interp >::get_fv,
	 "Return the y-axis array." )
    .def("__call__", py::vectorize(&utl::interpolator< utl::lin_interp >::operator()),
	 "Evaluate the interpolator at x (vectorised).\n"
	 "\nParameters\n----------\nx : float\n    Query point.\n"
	 "\nReturns\n-------\nfloat\n    Interpolated value.",
	 py::arg("x") )
    .def("integrate", &utl::interpolator< utl::lin_interp >::integrate,
	 INTEGRATE_DOC, py::arg("aa"), py::arg("bb") );

  py::class_< utl::interpolator< utl::log_interp > >( m, "log_interp",
    "Log-space piecewise-linear interpolator built from two equal-length arrays.\n"
    "Interpolation is performed in log10(x) vs log10(y) space.\n"
    INTERP_INIT_DOC )
    .def(py::init< const std::vector< double > &, const std::vector< double > & >())
    .def("get_x", &utl::interpolator< utl::log_interp >::get_xv,
	 "Return the x-axis array." )
    .def("get_y", &utl::interpolator< utl::log_interp >::get_fv,
	 "Return the y-axis array." )
    .def("__call__", py::vectorize(&utl::interpolator< utl::log_interp >::operator()),
	 "Evaluate the interpolator at x (vectorised).\n"
	 "\nParameters\n----------\nx : float\n    Query point.\n"
	 "\nReturns\n-------\nfloat\n    Interpolated value.",
	 py::arg("x") )
    .def("integrate", &utl::interpolator< utl::log_interp >::integrate,
	 INTEGRATE_DOC, py::arg("aa"), py::arg("bb") );

}
