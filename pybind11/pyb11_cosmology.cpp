// PyBind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
// External includes
#include <vector>
// Internal includes
#include <cosmological_model.h>

namespace py = pybind11;

PYBIND11_MODULE( cosmology, m ) {

  py::class_< scam::cosmo_model >( m, "model" )
    .def(py::init<
	 const double, const double, const double,
	 const double, const double, const double,
	 const double, const double, const double,
	 const double, const double, const double,
	 const std::size_t>(),
  	 py::arg("Om_M") = 0.3, py::arg("Om_b") = 0.045, py::arg("Om_L") = 0.7,
  	 py::arg("Om_n") = 0.0, py::arg("Om_r") = 0.0, py::arg("Om_K") = 0.0,
  	 py::arg("hh") = 0.7, py::arg("sigma8") = 0.8,
  	 py::arg("w_0") = 0.0, py::arg("w_a") = 0.0,
  	 py::arg("zmin") = 0.0, py::arg("zmax") = 100.,
  	 py::arg("thin") = 1000,
  	 "Constructor of the class cosmology.model" )
    .def("Hz", py::vectorize(&scam::cosmo_model::H_z),
	 "Hubble parameter at given redshift", py::arg("zz") )
    .def("dH", py::vectorize(&scam::cosmo_model::d_H),
	 "Hubble distance at given redshift", py::arg("zz") )
    .def("dC", py::vectorize(&scam::cosmo_model::d_C),
	 "Comoving distance at given redshift", py::arg("zz") )
    .def("ddC", py::vectorize(&scam::cosmo_model::dd_C),
	 "Derivative of the comoving distance at given redshift", py::arg("zz") )
    .def("dA", py::vectorize(&scam::cosmo_model::d_A),
	 "Angular diameter distance at given redshift", py::arg("zz") )
    .def("comoving_volume_unit", py::vectorize(&scam::cosmo_model::comoving_volume_unit),
	 "Comoving volume per redshift bin per unit steradian (dV/dzdOmega)",
	 py::arg("zz") )
    .def("comoving_volume", py::vectorize(&scam::cosmo_model::comoving_volume),
	 "Comoving volume at given redshift", py::arg("zz") )
    .def("cosmic_time", py::vectorize(&scam::cosmo_model::cosmic_time),
	 "Age of the Universe at given redshift (i.e. at the time photons were emitted)",
	 py::arg("zz") )
    .def("comoving_critical_density", py::vectorize(&scam::cosmo_model::rho_crit_comoving),
	 "Critical density at given redshift in units of [Msol Mpc^-3 h^2]",
	 py::arg("zz") )
    .def("critical_density", py::vectorize(&scam::cosmo_model::rho_crit),
	 "Critical density at given redshift in units of [Msol Mpc^-3]",
	 py::arg("zz") )
    .def("OmegaM", py::vectorize(&scam::cosmo_model::OmegaM),
	 "Matter density parameter at given redshift",
	 py::arg("zz") )
    .def("Omegab", py::vectorize(&scam::cosmo_model::Omegab),
	 "Baryonic matter density parameter at given redshift",
	 py::arg("zz") )
    .def("deltac", py::vectorize(&scam::cosmo_model::deltac),
	 "Critical overdensity threshold for virial collapse at given redshift",
	 py::arg("zz") )
    .def("D", py::vectorize(&scam::cosmo_model::DD),
	 "Amplitude of the growing solution of matter fluctuations at given redshift",
	 py::arg("zz") )
    .def("gz", py::vectorize(&scam::cosmo_model::gz),
	 "Growth factor at given redshift (i.e. (1+z)*D(z))",
	 py::arg("zz") )
    .def_readonly("zmin", &scam::cosmo_model::z_min, "Minimum redshift of internal grid")
    .def_readonly("zmax", &scam::cosmo_model::z_max, "Maximum redshift of internal grid");

} // end PYBIND11_MODULE
