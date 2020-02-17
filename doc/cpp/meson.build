#######################################################################
# External dependencies

# Doxygen executable

doxy_exe = find_program( [ 'doxygen' ],
                         required : true )

#######################################################################
# Doxygen documentation

man_path = '/'.join( [ get_option( 'prefix' ),
                       get_option( 'mandir' ) ] )

conf_data = configuration_data()
conf_data.set( 'DOXYGEN_OUTPUT_DIR', 'doc/cpp' )
conf_data.set( 'UTLTS_PATH', root + '/utiltities')
conf_data.set( 'INTRP_PATH', root + '/interpolator')
conf_data.set( 'COSMO_PATH', root + '/cosmology')
conf_data.set( 'FFTLW_PATH', root + '/fftlog_wrap')
conf_data.set( 'HMODL_PATH', root + '/halo_model')

doxy_file = configure_file( output : 'doxyfile',
                            input : 'doxyfile.in',
                            configuration : conf_data,
                            install : false )

doxy_target = custom_target( 'scampy_doxy',
                             input : doxy_file,
                             output : 'xml',
                             command : [ doxy_exe, doxy_file ],
                             install : true,
                             install_dir : man_path )

#######################################################################
