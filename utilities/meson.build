#######################################################################
# Internal includes

inc_dir = 'include'

inc_utl = include_directories( inc_dir )

#######################################################################
# Utilities Lib

utl_lib = library( 'utilities',
                   [ 'src/utilities.cpp',
                     'src/bit_manipulator.cpp' ],
       	 	   include_directories : inc_utl,
		   dependencies : [ gsl_dep ],
                   install : true  )

utl_dep = declare_dependency( link_with : utl_lib,
                              include_directories : inc_utl )

lib_headers += [ root + '/utilities/' + inc_dir + '/utilities.h',
                 root + '/utilities/' + inc_dir + '/error_handling.h',
                 root + '/utilities/' + inc_dir + '/bit_manipulator.h' ]

#######################################################################
