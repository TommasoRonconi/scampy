#######################################################################
# Internal includes

inc_dir = 'include'

inc_cosmo = include_directories( inc_dir )

#######################################################################
# Cosmology Lib

cosmo_lib = library( 'cosmo',
                     [ 'src/cosmology_interface.cpp',
                       'src/cosmology.cpp',
                       'src/luminosity_function.cpp' ],
       	 	     include_directories : [ inc_cosmo,
                                             inc_intrp ],
		     dependencies : [ gsl_dep,
                                      utl_dep,
                                      intrp_dep ],
                     install : true )

cosmo_dep = declare_dependency( link_with : cosmo_lib,
                                include_directories : [ inc_cosmo,
                                                        inc_intrp ] )

lib_headers += [ root + '/cosmology/' + inc_dir + '/cosmology_interface.h',
                 root + '/cosmology/' + inc_dir + '/cosmology.h' ]

#######################################################################
# Cosmology Lib C interface

cosmo_c_lib = library( 'cosmo_c_wrap',
                       [ 'src/cosmo_c_interface.cpp' ],
                       include_directories : [ inc_cosmo,
                                               inc_utl ],
                       dependencies : [ cosmo_dep,
                                        gsl_dep,
                                        utl_dep,
                                        intrp_dep ],
                       install : true )

cosmo_c_dep = declare_dependency( link_with : cosmo_c_lib,
                                  include_directories : [ inc_cosmo ] )

cosmo_c_path = custom_prefix + '/lib/libcosmo_c_wrap.so'

lib_headers += [ root + '/cosmology/' + inc_dir + '/cosmo_c_interface.h' ]

#######################################################################

