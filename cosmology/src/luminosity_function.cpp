#include <cosmology.h>
#include <gsl/gsl_sf_gamma.h>

//==============================================================================================

double scam::cosmology::phi_Bouwens15 ( const double LL, const double zz ) {

  const double phi_star = 0.47 * std::pow( 10, -0.27 * ( zz - 6. ) ) * 1.e-3;
  const double L_star = - 20.95 + 0.01 * ( zz - 6. );
  const double alpha = - 1.87 - 0.1 * ( zz - 6. ) + 1;
  const double expon = - 0.4 * ( LL - L_star );
  
  // Parametric Schechter Luminosity Function:
  return phi_star * 0.4 * scam::ln_10 * std::pow( 10, expon * alpha ) * std::exp( - std::pow( 10, expon ) );

}

//==============================================================================================

double scam::cosmology::phi_Bouwens16 ( const double LL, const double zz ) {

  const double phi_star = 0.45 * std::pow( 10, -0.21 * ( zz - 6. ) ) * 1.e-3;
  const double L_star = - 20.97 + 0.17 * ( zz - 6. );
  const double alpha = - 1.91 - 0.13 * ( zz - 6. ) + 1;
  const double expon = - 0.4 * ( LL - L_star );
  
  // Parametric Schechter Luminosity Function:
  return phi_star * 0.4 * scam::ln_10 * std::pow( 10, expon * alpha ) * std::exp( - std::pow( 10, expon ) );

}

//==============================================================================================

// double scam::cosmology::phi_Lapi17_uv ( const double & LL, const double & zz ) {

//   const vector< double > philog = { -1.96, -1.60, 4.22, -5.23 };
//   const vector< double > alphap = { 1.11, 2.85, -6.18, 4.44 };
//   const vector< double > sfrlog = { 0.01, 2.85, 0.43, -1.70 };

//   return phi_Lapi17 ( LL, zz, philog, sfrlog, alphap );

// }

// //==============================================================================================

// double scam::cosmology::phi_Lapi17_uvir ( const double & LL, const double & zz ) {

//   const vector< double > philog = { -2.13, -8.90, 18.07, -11.77 };
//   const vector< double > alphap = { 1.12, 3.73, -7.80, 5.15 };
//   const vector< double > sfrlog = { 0.72, 8.56, -10.08, 2.54 };

//   return phi_Lapi17 ( LL, zz, philog, sfrlog, alphap );

// }

double scam::cosmology::phi_Lapi17 ( const double & LL, const double & zz,
				     const std::vector< double > & param ) {

  const double csi = std::log10( 1 + zz );
  //const double Norm = scam::utl::poly_3 ( csi, { param.begin() + 0, param.begin() + 3 } );
  const double Norm = param[ 0 ] + param[ 1 ] * csi + param[ 2 ] * csi * csi + param[ 3 ] * csi * csi * csi;
  const double alpha = param[ 4 ] + param[ 5 ] * csi + param[ 6 ] * csi * csi + param[ 7 ] * csi * csi * csi;
  double sfrc = param[ 8 ] + param[ 9 ] * csi + param[ 10 ] * csi * csi + param[ 11 ] * csi * csi * csi;

  sfrc = std::pow( 10, sfrc );
  const double sfr = std::pow( 10, - 0.4 * ( LL - 4.83 ) - 9.8 ) / sfrc;
  // const double sfr = std::pow( 10, - 0.4 * LL - 7.4 ) / sfrc;
  
  return - 0.4 * Norm * std::pow( sfr, 1 - alpha ) * std::exp( - sfr );

}

//==============================================================================================

double scam::cosmology::phi_Lapi17_uv ( const double & LL, const double & zz ) {

  const std::vector< double > param = { -1.96, -1.60, 4.22, -5.23, //philog
					1.11, 2.85, -6.18, 4.44,   //alphap
					0.01, 2.85, 0.43, -1.70 }; //sfrlog

  return phi_Lapi17 ( LL, zz, param );

}

//==============================================================================================

double scam::cosmology::phi_Lapi17_uvir ( const double & LL, const double & zz ) {

  const std::vector< double > param = { -2.13, -8.90, 18.07, -11.77, //philog
					1.12, 3.73, -7.80, 5.15,     //alphap
					0.72, 8.56, -10.08, 2.54 };  //sfrlog
  
  return phi_Lapi17 ( LL, zz, param );

}

//==============================================================================================

// double scam::cosmology::phi_Lapi17 ( const double & LL, const double & zz,
// 					const std::vector< double > & philog,
// 					const std::vector< double > & sfrlog,
// 					const std::vector< double > & alphap ) {

//   const double csi = std::log10( 1 + zz );
//   const double Norm = philog[ 0 ] + philog[ 1 ] * csi + philog[ 2 ] * csi * csi;
//   const double alpha = alphap[ 0 ] + alphap[ 1 ] * csi + alphap[ 2 ] * csi * csi;
//   double sfrc = sfrlog[ 0 ] + sfrlog[ 1 ] * csi + sfrlog[ 2 ] * csi * csi;

//   sfrc = std::pow( 10, sfrc );
//   const double sfr = std::pow( 10, - 0.4 * LL - 7.4 ) / sfrc;
  
//   return - 0.4 * Norm * std::pow( sfr, 1. - alpha ) * std::exp( - sfr );

// }

//==============================================================================================

double scam::cosmology::phi (const double & LL, const double & zz, const std::string & modelLF) {

  // Parameterisation valid in redshift range 4 <~ z <~ 8 (Bouwens et al., 2015):
  if ( modelLF == "Bouwens15")
    return phi_Bouwens15( LL, zz );

  // Parameterization valid in redshift range 4 <~ z <~ 10 (Bouwens et al., 2016):
  else if (modelLF == "Bouwens16")
    return phi_Bouwens16( LL, zz );
  
  // Parameterization from Lapi et al., 2016  
  else if (modelLF == "Lapi17_uv")
    return phi_Lapi17_uv( LL, zz );
  
  // Parameterization from Lapi et al., 2016  
  else if (modelLF == "Lapi17_uvir")
    return phi_Lapi17_uvir( LL, zz );
  
  else
    throw scam_err::type_invalid {
      "Parameterisation " + modelLF + " not accepted!\n" +
	"Accepted parameterisations are 'Bouwens15', 'Bouwens16', 'Lapi17_uv', 'Lapi17_uvir'."
	};
  
}

//==============================================================================================
