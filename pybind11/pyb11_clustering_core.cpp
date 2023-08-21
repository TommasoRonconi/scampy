// PyBind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// External includes
#include <vector>
// Internal includes
#include <clustering_core.h>

namespace py = pybind11;

// template< &utl::dist, const float d, const float box >
// std::vector< std::size_t > utl::d2D_DD_generic ();

PYBIND11_MODULE( clustering_core, m ) {

  // 2D block
  m.def( "d2D_DD", &utl::d2D_DD,
	 "Compute all distances among 2D coordinates",
	 py::arg("X"), py::arg("Y"), py::arg("rbin") );
  m.def( "d2D_DD_omp", &utl::d2D_DD_omp,
	 "Compute all distances among 2D coordinates",
	 py::arg("X"), py::arg("Y"), py::arg("rbin") );
  m.def( "d2D_DR", &utl::d2D_DR,
	 "Compute all distances among 2D coordinates "
	 "from two sets of coordinates",
	 py::arg("X1"), py::arg("Y1"),
	 py::arg("X2"), py::arg("Y2"),
	 py::arg("rbin") );
  m.def( "d2D_DR_omp", &utl::d2D_DR_omp,
	 "Compute all distances among 2D coordinates "
	 "from two sets of coordinates",
	 py::arg("X1"), py::arg("Y1"),
	 py::arg("X2"), py::arg("Y2"),
	 py::arg("rbin") );

  // 2D-Angular block
  m.def( "dA2D_DD", &utl::dA2D_DD,
	 "Compute all distances among 2D angular coordinates",
	 py::arg("RA"), py::arg("Dec"), py::arg("thetabin") );
  m.def( "dA2D_DD_omp", &utl::dA2D_DD_omp,
	 "Compute all distances among 2D angular coordinates",
	 py::arg("RA"), py::arg("Dec"), py::arg("thetabin") );
  m.def( "dA2D_DR", &utl::dA2D_DR,
	 "Compute all distances among 2D angular coordinates "
	 "from two sets of coordinates",
	 py::arg("RA1"), py::arg("Dec1"),
	 py::arg("RA2"), py::arg("Dec2"),
	 py::arg("thetabin") );
  m.def( "dA2D_DR_omp", &utl::dA2D_DR_omp,
	 "Compute all distances among 2D angular coordinates "
	 "from two sets of coordinates",
	 py::arg("RA1"), py::arg("Dec1"),
	 py::arg("RA2"), py::arg("Dec2"),
	 py::arg("thetabin") );

  // 3D block
  m.def( "d3D_DD", &utl::d3D_DD,
	 "Compute all distances among 3D coordinates",
	 py::arg("X"), py::arg("Y"), py::arg("Z"), py::arg("rbin") );
  m.def( "d3D_DD_omp", &utl::d3D_DD_omp,
	 "Compute all distances among 3D coordinates",
	 py::arg("X"), py::arg("Y"), py::arg("Z"), py::arg("rbin") );
  m.def( "d3D_DR", &utl::d3D_DR,
	 "Compute all distances among 3D coordinates "
	 "from two sets of coordinates",
	 py::arg("X1"), py::arg("Y1"), py::arg("Z1"),
	 py::arg("X2"), py::arg("Y2"), py::arg("Z2"),
	 py::arg("rbin") );
  m.def( "d3D_DR_omp", &utl::d3D_DR_omp,
	 "Compute all distances among 3D coordinates "
	 "from two sets of coordinates",
	 py::arg("X1"), py::arg("Y1"), py::arg("Z1"),
	 py::arg("X2"), py::arg("Y2"), py::arg("Z2"),
	 py::arg("rbin") );

}
