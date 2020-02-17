#######################################################################
# External dependencies

subdir( 'cpp' )

# Sphynx executable

sphinx_exe = find_program( [ 'sphinx-build' ],
                           required : true )

#######################################################################
# Sphinx documentation

sphinx_target = custom_target( 'scampy_sphinx',
                               output : 'html',
                               command : [ sphinx_exe,
                                           '-b', 'html',
                                           '-Dbreathe_projects.ScamPy=' +
                                           meson.build_root() +
                                           '/doc/cpp/xml',
                                           meson.current_source_dir(),
                                           'doc/html' ],
                               install : true,
                               install_dir : man_path )

#######################################################################
