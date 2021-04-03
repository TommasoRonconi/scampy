#######################################################################
# Search for python installation

pymod = import( 'python' )
py = pymod.find_installation()
message( py.get_install_dir() )

#######################################################################
# Configure python module

conf_data = configuration_data()
conf_data.set( 'INTRP_WRAP_PATH', intrp_c_path )
conf_data.set( 'COSMO_WRAP_PATH', cosmo_c_path )
conf_data.set( 'OCP_WRAP_PATH', ocp_c_path )
conf_data.set( 'HM_WRAP_PATH', hm_c_path )

configure_file( output : 'cwrap.py',
                input : 'cwrap.py.in',
                configuration : conf_data,
                install : true,
                install_dir :  py.get_install_dir() + 'scampy/cwrap' )

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

py.install_sources( [ 'scampy/cwrap/__init__.py' ],
                    subdir : 'scampy/cwrap' )

#######################################################################


