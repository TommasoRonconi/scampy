#######################################################################
# C++ tests

tests_cxx = [ 'occupation_p', 'halo_model' ]

foreach _t : tests_cxx

exe_cxx = executable( 'test_' + _t + '.x',
                      [ 'test_' + _t + '.cpp' ],
                      include_directories : [ inc_hm,
                                              inc_utl ],
                      dependencies : [ hm_dep,
                                       cosmo_dep,
                                       gsl_dep,
                                       utl_dep ] )

# test( _t + '_run_and_diff',
#       compare,
#       args : [ exe.full_path(),
#                _prefix + '/output/' + _t + '_ref.out' ],
#       workdir : meson.current_build_dir() )

endforeach

#######################################################################
# C tests

# tests_c = [ 'cosmo_c_wrap', 'hm_c_wrap' ]

# foreach _t : tests_c

# exe_c = executable( 'test_' + _t + '.x',
#                     [ 'test_' + _t + '.c' ],
#                     include_directories : [ inc_hm,
#                                             inc_utl ],
#                     dependencies : [ cosmo_c_dep,
#                                      hm_c_dep,
#                                      ocp_c_dep,
#                                      gsl_dep,
#                                      utl_dep ],
#                     install : true,
#                     override_options : [ 'compiler=c', 'c_std=c99' ] )

# # test( _t + '_run_and_diff',
# #       compare,
# #       args : [ exe.full_path(),
# #                _prefix + '/output/' + _t + '_ref.out' ],
# #       workdir : meson.current_build_dir() )

# endforeach

#######################################################################
