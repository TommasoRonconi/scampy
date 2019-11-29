#######################################################################
# Google unit testing

gtest = dependency('gtest',
                   fallback : ['gtest', 'gtest_main_dep'],  # magic of meson wrap
                   required: true)

_tests = [ 'occupation_p' ]
# _tests =  [ 'mathematics_01',
#        	    'physics_01' ]

foreach _t : _tests
	_e = executable( _t, _t+'.cpp',
	     		 dependencies: [ hm_dep,
                                         cosmo_dep,
                                         gsl_dep,
                                         utl_dep,
                                         gtest ] )
	test( _t, _e )

endforeach

#######################################################################
