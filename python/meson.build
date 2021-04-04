#######################################################################
# Configure python module

conf_data = configuration_data()

if get_option( 'full-build' )
infile = 'cwrap-fb1.py.in'
conf_data.set( 'INTRP_WRAP_PATH', intrp_c_path )
conf_data.set( 'COSMO_WRAP_PATH', cosmo_c_path )
conf_data.set( 'OCP_WRAP_PATH', ocp_c_path )
conf_data.set( 'HM_WRAP_PATH', hm_c_path )
wrapdir = '/'.join( [ py.get_install_dir(), 'scampy', 'internal' ] )
else
infile = 'cwrap-fb0.py.in'
conf_data.set( 'INTRP_WRAP_PATH', '/' + intrp_c_lib_name )
conf_data.set( 'COSMO_WRAP_PATH', '/' + cosmo_c_lib_name )
conf_data.set( 'OCP_WRAP_PATH', '/' + ocp_c_lib_name )
conf_data.set( 'HM_WRAP_PATH', '/' + hm_c_lib_name )
wrapdir = libdir
endif

configure_file( output : 'cwrap.py',
                input : infile,
                configuration : conf_data,
                install : true,
                install_dir : wrapdir )

py.install_sources( [ 'scampy/__init__.py',
                      'scampy/abundance_matching.py',
                      'scampy/catalogue.py',
                      'scampy/cosmology.py',
                      'scampy/gadget_file.py',
                      'scampy/halo_model.py',
                      'scampy/interpolator.py',
                      'scampy/objects.py',
                      'scampy/occupation_p.py' ],
                    subdir : 'scampy' )

py.install_sources( [ 'scampy/internal/__init__.py' ],
                    subdir : 'scampy/internal' )

#######################################################################


