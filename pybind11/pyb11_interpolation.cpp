#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <interpolation.h>
#include <vector>

namespace py = pybind11;

template class utl::interpolator< utl::lin_interp >;

PYBIND11_MODULE( interpolation, m ) {
  
  py::class_< utl::interpolator< utl::lin_interp > >( m, "lin_interp" )
    .def(py::init< const std::vector< double > &, const std::vector< double > & >())
    .def("get_x", &utl::interpolator< utl::lin_interp >::get_xv )
    .def("get_y", &utl::interpolator< utl::lin_interp >::get_fv )
    .def("__call__", py::vectorize(&utl::interpolator< utl::lin_interp >::operator()),
	 "Evaluate the interpolated function", py::arg("x") )
    .def("integrate", &utl::interpolator< utl::lin_interp >::integrate,
	 "Integrate the interpolated function within a specified interval.\n"
	 "\nParameters"
	 "\n----------"
	 "\naa : float"
	 "\n\tinterval lower limit"
	 "\nbb : float"
	 "\n\tinterval upper limit\n"
	 "\nReturns"
	 "\n------"
	 "\n: float"
	 "\n\tapproximate integral", py::arg("aa"), py::arg("bb") );
  
}
