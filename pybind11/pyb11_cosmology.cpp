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
  	 py::arg("w_0") = -1.0, py::arg("w_a") = 0.0,
  	 py::arg("zmin") = 0.0, py::arg("zmax") = 100.,
  	 py::arg("thin") = 1000,
	 "Flat or curved FLRW cosmological model.\n"
	 "\nParameters\n"
	 "----------\n"
	 "Om_M : float, optional\n"
	 "    Total matter density parameter Omega_M (default: 0.3).\n"
	 "Om_b : float, optional\n"
	 "    Baryonic matter density parameter Omega_b (default: 0.045).\n"
	 "Om_L : float, optional\n"
	 "    Dark-energy density parameter Omega_Lambda (default: 0.7).\n"
	 "Om_n : float, optional\n"
	 "    Massive-neutrino density parameter (default: 0.0).\n"
	 "Om_r : float, optional\n"
	 "    Radiation density parameter (default: 0.0).\n"
	 "Om_K : float, optional\n"
	 "    Curvature density parameter (default: 0.0).\n"
	 "hh : float, optional\n"
	 "    Dimensionless Hubble constant h = H0/100 (default: 0.7).\n"
	 "sigma8 : float, optional\n"
	 "    RMS matter-fluctuation amplitude at 8 Mpc/h (default: 0.8).\n"
	 "w_0 : float, optional\n"
	 "    Dark-energy equation-of-state constant term (default: -1.0).\n"
	 "w_a : float, optional\n"
	 "    Dark-energy equation-of-state linear term w(a) = w_0 + w_a*(1-a) (default: 0.0).\n"
	 "zmin : float, optional\n"
	 "    Minimum redshift of the internal tabulation grid (default: 0.0).\n"
	 "zmax : float, optional\n"
	 "    Maximum redshift of the internal tabulation grid (default: 100.0).\n"
	 "thin : int, optional\n"
	 "    Number of points in the internal tabulation grid (default: 1000)." )
    .def("Hz", py::vectorize(&scam::cosmo_model::H_z),
	 "Hubble parameter H(z) in km/s/Mpc.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    H(z) [km/s/Mpc].", py::arg("zz") )
    .def("dH", py::vectorize(&scam::cosmo_model::d_H),
	 "Hubble distance c/H(z) in Mpc.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    c/H(z) [Mpc].", py::arg("zz") )
    .def("dC", py::vectorize(&scam::cosmo_model::d_C),
	 "Line-of-sight comoving distance in Mpc.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    d_C(z) [Mpc].", py::arg("zz") )
    .def("ddC", py::vectorize(&scam::cosmo_model::dd_C),
	 "Derivative of the comoving distance with respect to redshift, dd_C/dz.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    dd_C/dz [Mpc].", py::arg("zz") )
    .def("dA", py::vectorize(&scam::cosmo_model::d_A),
	 "Angular diameter distance in Mpc.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    d_A(z) [Mpc].", py::arg("zz") )
    .def("comoving_volume_unit", py::vectorize(&scam::cosmo_model::comoving_volume_unit),
	 "Comoving volume element per unit redshift per unit steradian, dV/(dz dOmega).\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    dV/(dz dOmega) [(Mpc/h)^3 / sr].",
	 py::arg("zz") )
    .def("comoving_volume", py::vectorize(&scam::cosmo_model::comoving_volume),
	 "Total comoving volume enclosed within redshift z in (Mpc/h)^3.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    V_C(z) [(Mpc/h)^3].", py::arg("zz") )
    .def("cosmic_time", py::vectorize(&scam::cosmo_model::cosmic_time),
	 "Age of the Universe at redshift z (lookback time from z to infinity) in yr.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    t(z) [yr].",
	 py::arg("zz") )
    .def("critical_density_comoving", py::vectorize(&scam::cosmo_model::rho_crit_comoving),
	 "Critical density in comoving units [Msol Mpc^-3 h^2].\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    rho_crit(z) [Msol Mpc^-3 h^2].",
	 py::arg("zz") )
    .def("critical_density", py::vectorize(&scam::cosmo_model::rho_crit),
	 "Critical density in physical units [Msol Mpc^-3].\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    rho_crit(z) [Msol Mpc^-3].",
	 py::arg("zz") )
    .def("OmegaM", py::vectorize(&scam::cosmo_model::OmegaM),
	 "Total matter density parameter Omega_M(z) = rho_M(z)/rho_crit(z).\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    Omega_M(z).", py::arg("zz") )
    .def("Omegab", py::vectorize(&scam::cosmo_model::Omegab),
	 "Baryonic matter density parameter Omega_b(z).\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    Omega_b(z).", py::arg("zz") )
    .def("deltac", py::vectorize(&scam::cosmo_model::deltac),
	 "Linear critical overdensity for spherical collapse delta_c(z).\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    delta_c(z).", py::arg("zz") )
    .def("D", py::vectorize(&scam::cosmo_model::DD),
	 "Linear growth factor D(z), normalised to D(0) = 1.\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    D(z).", py::arg("zz") )
    .def("gz", py::vectorize(&scam::cosmo_model::gz),
	 "Unnormalised growth factor g(z) = (1+z)*D(z).\n"
	 "\nParameters\n----------\nzz : float\n    Redshift.\n"
	 "\nReturns\n-------\nfloat\n    (1+z)*D(z).", py::arg("zz") )
    .def("Delta_c", py::vectorize(&scam::cosmo_model::Delta_c_NakamuraSuto98),
	 "Virial overdensity Delta_c(z) from Nakamura & Suto (1998).\n"
	 "\nReturns\n-------\nfloat\n    Delta_c(z)." )
    .def_readonly("H0", &scam::cosmo_model::H0,
		  "Hubble constant at z=0 (defined as 100*h0 and expressed in km/s/Mpc)")
    .def_readonly("t_H0", &scam::cosmo_model::t_H0,
		  "Hubble time (defined as 1/H0 and expressed in yr)")
    .def_readonly("d_H0", &scam::cosmo_model::d_H0,
		  "Hubble orizon distance (defined as 1/H0 and expressed in Mpc)")
    .def_readonly("zmin", &scam::cosmo_model::z_min, "Minimum redshift of internal grid")
    .def_readonly("zmax", &scam::cosmo_model::z_max, "Maximum redshift of internal grid")
    .def_readonly("param", &scam::cosmo_model::param, "Cosmological parameters");

} // end PYBIND11_MODULE
