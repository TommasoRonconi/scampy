#######################################################################
# Google unit testing

gtest = dependency('gtest',
                   fallback : ['gtest', 'gtest_main_dep'],
                   required: true)

_tests = [ 'test_interpolator',
           'test_occupation_p',
           'test_cosmology' ]

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
