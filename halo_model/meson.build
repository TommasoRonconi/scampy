#######################################################################
# Internal includes

inc_dir = 'include'

inc_hm = include_directories( inc_dir )

#######################################################################
# Halo-Model Lib

if get_option( 'full-build' )
hm_lib = library( 'halo_model',
                  [ 'src/halo_model.cpp',
                    'src/cross_halo_model.cpp',
                    'src/occupation_p.cpp',
                    'src/harikane16_p.cpp',
                    'src/tinker10_p.cpp' ],
       	 	  include_directories : [ inc_hm,
                                          inc_cosmo,
                                          inc_intrp,
                                          inc_fftlwrap,
                                          inc_utl ],
		  dependencies : [ omp_dep,
                                   gsl_dep,
                                   utl_dep,
                                   intrp_dep,
                                   fftlwrap_dep,
                                   cosmo_dep ],
                  install : true,
                  install_dir : libdir )
else
hm_lib = static_library( 'halo_model',
                         [ 'src/halo_model.cpp',
                           'src/cross_halo_model.cpp',
                           'src/occupation_p.cpp',
                           'src/harikane16_p.cpp',
                           'src/tinker10_p.cpp' ],
       	 	         include_directories : [ inc_hm,
                                                 inc_cosmo,
                                                 inc_intrp,
                                                 inc_fftlwrap,
                                                 inc_utl ],
		         dependencies : [ omp_dep,
                                          gsl_dep,
                                          utl_dep,
                                          intrp_dep,
                                          fftlwrap_dep,
                                          cosmo_dep ],
                         install : true,
                         install_dir : libdir )
endif

hm_dep = declare_dependency( link_with : hm_lib,
                             include_directories : [ inc_hm,
                                                     inc_cosmo,
                                                     inc_intrp,
                                                     inc_fftlwrap,
                                                     inc_utl] )

#######################################################################
# Occupation C Interface

ocp_c_lib = library( 'ocp_c_wrap',
                     [ 'src/occupation_c_interface.cpp' ],
                     include_directories : [ inc_hm,
                                             inc_cosmo,
                                             inc_intrp,
                                             inc_fftlwrap ],
		     dependencies : [ hm_dep,
                                      omp_dep,
                                      gsl_dep,
                                      utl_dep,
                                      intrp_dep,
                                      fftlwrap_dep,
                                      cosmo_dep ],
                     install : true,
                     install_dir : libdir )

ocp_c_dep = declare_dependency( link_with : ocp_c_lib,
                                include_directories : [ inc_hm ] )

ocp_c_lib_name = ocp_c_lib.full_path().split('/')[-1]

ocp_c_path = '/'.join( [ get_option( 'prefix' ),
                         get_option( 'libdir' ),
                         ocp_c_lib_name ] )

#######################################################################
# Halo-Model Lib C Interface

hm_c_lib = library( 'hm_c_wrap',
                    [ 'src/halo_model_c_interface.cpp' ],
                    include_directories : [ inc_hm,
                                            inc_cosmo,
                                            inc_intrp,
                                            inc_fftlwrap ],
		    dependencies : [ hm_dep,
                                     ocp_c_dep,
                                     omp_dep,
                                     gsl_dep,
                                     utl_dep,
                                     intrp_dep,
                                     fftlwrap_dep,
                                     cosmo_dep ],
                    install : true,
                    install_dir : libdir )

hm_c_dep = declare_dependency( link_with : hm_c_lib,
                               include_directories : [ inc_hm ] )

hm_c_lib_name = hm_c_lib.full_path().split('/')[-1]

hm_c_path = '/'.join( [ get_option( 'prefix' ),
                        get_option( 'libdir' ),
                        hm_c_lib_name ] )

#######################################################################
# Add headers to list

if get_option( 'full-build' )
lib_headers += [ root + '/halo_model/' + inc_dir + '/halo_model.h',
                 root + '/halo_model/' + inc_dir + '/cross_halo_model.h',
                 root + '/halo_model/' + inc_dir + '/occupation_p.h',
                 root + '/halo_model/' + inc_dir + '/harikane16_p.h',
                 root + '/halo_model/' + inc_dir + '/tinker10_p.h',
                 root + '/halo_model/' + inc_dir + '/occupation_c_interface.h',
                 root + '/halo_model/' + inc_dir + '/halo_model_c_interface.h' ]
endif

#######################################################################
