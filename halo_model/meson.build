#######################################################################
# Internal includes

inc_dir = 'include'

inc_hm = include_directories( inc_dir )

#######################################################################
# Halo-Model Lib

hm_lib = library( 'halo_model',
                  [ 'src/halo_model.cpp',
                    'src/occupation_p.cpp',
                    'src/harikane16_p.cpp',
                    'src/tinker10_p.cpp' ],
       	 	  include_directories : [ inc_hm,
                                          inc_cosmo,
                                          inc_intrp,
                                          inc_fftlwrap ],
		  dependencies : [ omp_dep,
                                   gsl_dep,
                                   utl_dep,
                                   intrp_dep,
                                   fftlwrap_dep,
                                   cosmo_dep ],
                  install : true )

hm_dep = declare_dependency( link_with : hm_lib,
                             include_directories : [ inc_hm,
                                                     inc_cosmo,
                                                     inc_intrp,
                                                     inc_fftlwrap ] )

lib_headers += [ root + '/halo_model/' + inc_dir + '/halo_model.h' ]
lib_headers += [ root + '/halo_model/' + inc_dir + '/occupation_p.h' ]
lib_headers += [ root + '/halo_model/' + inc_dir + '/harikane16_p.h' ]
lib_headers += [ root + '/halo_model/' + inc_dir + '/tinker10_p.h' ]

#######################################################################
# Halo-Model C Lib Interface

hm_c_lib = library( 'hm_c_wrap',
                    [ 'src/halo_model_c_interface.cpp' ],
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
                    install : true )

hm_c_dep = declare_dependency( link_with : hm_c_lib,
                               include_directories : [ inc_hm ] )

hm_c_lib_name = hm_c_lib.full_path().split('/')[-1]

hm_c_path = '/'.join( [ get_option( 'prefix' ),
                        get_option( 'libdir' ),
                        hm_c_lib_name ] )

lib_headers += [ root + '/halo_model/' + inc_dir + '/halo_model_c_interface.h' ]

#######################################################################
