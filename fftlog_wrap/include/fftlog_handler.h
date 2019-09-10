/**
 *  @file include/fftlog_handler.h
 *
 *  @brief Header to handle with FFTLog wrapping  
 *
 *  In this file custom types and functions to handle
 *  Fast Fourier Transforms in log-space via the FFTLog 
 *  library are declared.
 *
 *  @author Tommaso Ronconi
 *
 *  @author tronconi@sissa.it
 */


#ifndef __FFTLOG_HANDLER__
#define __FFTLOG_HANDLER__

#include <utilities.h>
#include <interpolation.h>

namespace sico {

  extern "C" {

    /**
     *  @brief wrapper of the fhti subroutine contained in
     *  External/fftlog-f90-master/fftlog.f This is an
     *  initialization routine
     *
     *  @param [in] _n number of points in the array to be transformed
     *  @param [in] _mu index of \f$J_\mu\f$ in Hankel transform
     *  @param [in] _q exponent of power law bias
     *  @param [in] _dlnr separation between natural log of points
     *
     *  @param [in] _kr \f$k_c\cdot r_c\f$ where c is central point
     *  of array \f$kr = k_j r_(n+1-j) = k_(n+1-j) r_j\f$
     *
     *  @param [in] _kropt 0 &rarr; use input kr as is; 1 &rarr;
     *  change kr to nearest low-ringing kr, quietly; 2 &rarr;
     *  change kr to nearest low-ringing kr, verbosely; 3 &rarr;
     *  change kr interactively
     *
     *  @param [out] _wsave working array
     *  @param [out] _ok 1 &rarr; all went ok; 0 &rarr; error in initializations
     *
     *  @return none
     */
    void fhti_ ( int * size, double * mu,
		 double * qq, double * dlnr,
		 double * kr, int * kropt,
		 double * wsave, int * ok);
    
    void fht_ ( int * size, double * ap,
		int * dir, double * wsave);

  }

  namespace utl {

    class fftlog_handler {

    private:

      int _kropt = 1;
      int _ok;
      double _dlnx;
      double * _wsave;

    protected:

      int _dir;
      int _size;
      double _kr = 1.;
      double _q = 0.;
      double _mu = 0.;
      double _ci, _lnx_med, _lnk_med;
      double * _ap;
      std::vector< double > _xx, _fx;
      std::vector< double > _kk, _fk;
      void _transform () { fht_( &_size, _ap, &_dir, _wsave ); }
      
    public:
      
      interpolator< gsl_log_interp > fk;

      fftlog_handler () {

      	_ap = new double[ 128 ];
	_wsave = new double[ 512 ];

      }

      fftlog_handler ( const std::vector< double> & xx,
		       const std::vector< double> & fx,
		       const double kr = 1., const int dir = 1,
		       const double qq = 0., const double mu = 0.5 );

      virtual ~fftlog_handler () {

	delete[] _wsave;
      	delete[] _ap;
	
      }

      virtual void transform () = 0;

      virtual std::vector< double > transform ( const std::vector< double > & kk ) = 0;

    }; // endclass fftlog_handler

    class fftlog_3Dspace : public fftlog_handler {

    private:

      const double _fact = 0.5 * sico::ip * sico::ip * std::sqrt( 0.5 * sico::pi );

    public:

      fftlog_3Dspace ( const std::vector< double> & xx,
		       const std::vector< double> & fx,
		       const double kr = 1. );

      ~fftlog_3Dspace () = default;

      void transform () override;

      std::vector< double > transform ( const std::vector< double > & kk ) override;

    }; // endclass fftlog_3Dspace

    class fftlog_projected : public fftlog_handler  {

    private:

      const double _fact = 0.5 * sico::ip; 

    public:

      fftlog_projected ( const std::vector< double> & xx,
			 const std::vector< double> & fx,
			 const double kr = 1. );

      ~fftlog_projected () = default;

      void transform () override;

      std::vector< double > transform ( const std::vector< double > & kk ) override;
      
    }; // endclass fftlog_projected
    
  } // endnamespace utl

} // endnamespace sico

#endif //__FFTLOG_HANDLER__
