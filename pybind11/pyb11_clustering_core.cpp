// PyBind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// External includes
#include <vector>
// Internal includes
#include <clustering_core.h>

namespace py = pybind11;

#define DD2D_DOC \
  "Count data-data pairs in 2D separation bins.\n" \
  "\nParameters\n----------\n" \
  "X : list of float\n    X-coordinates of the catalogue.\n" \
  "Y : list of float\n    Y-coordinates of the catalogue.\n" \
  "rbin : list of float\n    Separation bin edges.\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

#define DR2D_DOC \
  "Count data-random cross-pairs in 2D separation bins.\n" \
  "\nParameters\n----------\n" \
  "X1 : list of float\n    X-coordinates of the first (data) catalogue.\n" \
  "Y1 : list of float\n    Y-coordinates of the first catalogue.\n" \
  "X2 : list of float\n    X-coordinates of the second (random) catalogue.\n" \
  "Y2 : list of float\n    Y-coordinates of the second catalogue.\n" \
  "rbin : list of float\n    Separation bin edges.\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

#define DA2D_DD_DOC \
  "Count data-data pairs in 2D angular separation bins.\n" \
  "\nParameters\n----------\n" \
  "RA : list of float\n    Right-ascension coordinates [rad].\n" \
  "Dec : list of float\n    Declination coordinates [rad].\n" \
  "thetabin : list of float\n    Angular separation bin edges [rad].\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

#define DA2D_DR_DOC \
  "Count data-random cross-pairs in 2D angular separation bins.\n" \
  "\nParameters\n----------\n" \
  "RA1 : list of float\n    Right-ascension of the first (data) catalogue [rad].\n" \
  "Dec1 : list of float\n    Declination of the first catalogue [rad].\n" \
  "RA2 : list of float\n    Right-ascension of the second (random) catalogue [rad].\n" \
  "Dec2 : list of float\n    Declination of the second catalogue [rad].\n" \
  "thetabin : list of float\n    Angular separation bin edges [rad].\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

#define DD3D_DOC \
  "Count data-data pairs in 3D separation bins.\n" \
  "\nParameters\n----------\n" \
  "X : list of float\n    X-coordinates of the catalogue.\n" \
  "Y : list of float\n    Y-coordinates of the catalogue.\n" \
  "Z : list of float\n    Z-coordinates of the catalogue.\n" \
  "rbin : list of float\n    Separation bin edges.\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

#define DR3D_DOC \
  "Count data-random cross-pairs in 3D separation bins.\n" \
  "\nParameters\n----------\n" \
  "X1 : list of float\n    X-coordinates of the first (data) catalogue.\n" \
  "Y1 : list of float\n    Y-coordinates of the first catalogue.\n" \
  "Z1 : list of float\n    Z-coordinates of the first catalogue.\n" \
  "X2 : list of float\n    X-coordinates of the second (random) catalogue.\n" \
  "Y2 : list of float\n    Y-coordinates of the second catalogue.\n" \
  "Z2 : list of float\n    Z-coordinates of the second catalogue.\n" \
  "rbin : list of float\n    Separation bin edges.\n" \
  "\nReturns\n-------\nlist of int\n    Pair counts per bin."

PYBIND11_MODULE( clustering_core, m ) {

  // 2D block
  m.def( "d2D_DD", &utl::d2D_DD, DD2D_DOC,
	 py::arg("X"), py::arg("Y"), py::arg("rbin") );
  m.def( "d2D_DD_omp", &utl::d2D_DD_omp,
	 DD2D_DOC " Uses OpenMP parallelism.",
	 py::arg("X"), py::arg("Y"), py::arg("rbin") );
  m.def( "d2D_DR", &utl::d2D_DR, DR2D_DOC,
	 py::arg("X1"), py::arg("Y1"),
	 py::arg("X2"), py::arg("Y2"),
	 py::arg("rbin") );
  m.def( "d2D_DR_omp", &utl::d2D_DR_omp,
	 DR2D_DOC " Uses OpenMP parallelism.",
	 py::arg("X1"), py::arg("Y1"),
	 py::arg("X2"), py::arg("Y2"),
	 py::arg("rbin") );

  // 2D-Angular block
  m.def( "dA2D_DD", &utl::dA2D_DD, DA2D_DD_DOC,
	 py::arg("RA"), py::arg("Dec"), py::arg("thetabin") );
  m.def( "dA2D_DD_omp", &utl::dA2D_DD_omp,
	 DA2D_DD_DOC " Uses OpenMP parallelism.",
	 py::arg("RA"), py::arg("Dec"), py::arg("thetabin") );
  m.def( "dA2D_DR", &utl::dA2D_DR, DA2D_DR_DOC,
	 py::arg("RA1"), py::arg("Dec1"),
	 py::arg("RA2"), py::arg("Dec2"),
	 py::arg("thetabin") );
  m.def( "dA2D_DR_omp", &utl::dA2D_DR_omp,
	 DA2D_DR_DOC " Uses OpenMP parallelism.",
	 py::arg("RA1"), py::arg("Dec1"),
	 py::arg("RA2"), py::arg("Dec2"),
	 py::arg("thetabin") );

  // 3D block
  m.def( "d3D_DD", &utl::d3D_DD, DD3D_DOC,
	 py::arg("X"), py::arg("Y"), py::arg("Z"), py::arg("rbin") );
  m.def( "d3D_DD_omp", &utl::d3D_DD_omp,
	 DD3D_DOC " Uses OpenMP parallelism.",
	 py::arg("X"), py::arg("Y"), py::arg("Z"), py::arg("rbin") );
  m.def( "d3D_DR", &utl::d3D_DR, DR3D_DOC,
	 py::arg("X1"), py::arg("Y1"), py::arg("Z1"),
	 py::arg("X2"), py::arg("Y2"), py::arg("Z2"),
	 py::arg("rbin") );
  m.def( "d3D_DR_omp", &utl::d3D_DR_omp,
	 DR3D_DOC " Uses OpenMP parallelism.",
	 py::arg("X1"), py::arg("Y1"), py::arg("Z1"),
	 py::arg("X2"), py::arg("Y2"), py::arg("Z2"),
	 py::arg("rbin") );

}
