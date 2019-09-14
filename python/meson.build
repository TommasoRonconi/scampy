#######################################################################
# Install python module

install_subdir( 'mockpy', install_dir : 'python' )

#######################################################################
# Configure python module

conf_data = configuration_data()
conf_data.set( 'INTRP_WRAP_PATH', intrp_c_path )
conf_data.set( 'COSMO_WRAP_PATH', cosmo_c_path )
conf_data.set( 'HM_WRAP_PATH', hm_c_path )

configure_file( output : 'interpolator.py',
                input : 'interpolator.py.in',
                configuration : conf_data,
                install : true,
                install_dir : 'python/mockpy' )

configure_file( output : 'cosmology.py',
                input : 'cosmology.py.in',
                configuration : conf_data,
                install : true,
                install_dir : 'python/mockpy' )

configure_file( output : 'halo_model.py',
                input : 'halo_model.py.in',
                configuration : conf_data,
                install : true,
                install_dir : 'python/mockpy' )

#######################################################################


