tests = [ 'occupation_p' ]

foreach _t : tests

exe = executable( 'test_' + _t + '.x',
                  [ 'test_' + _t + '.cpp' ],
                  include_directories : [ inc_hm,
                                          inc_utl ],
                  dependencies : [ hm_dep,
                                   gsl_dep,
                                   utl_dep ],
                  install : true )

# test( _t + '_run_and_diff',
#       compare,
#       args : [ exe.full_path(),
#                _prefix + '/output/' + _t + '_ref.out' ],
#       workdir : meson.current_build_dir() )

endforeach
