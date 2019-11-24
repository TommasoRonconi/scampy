#######################################################################
# Sphinx documentation

sphinx_target = custom_target( 'scampy_sphinx',
                               output : 'html',
                               command : [ sphinx_exe,
                                           '-b', 'html',
                                           '-Dbreathe_projects.ScamPy=' +
                                           meson.build_root() +
                                           '/doc/xml',
                                           meson.current_source_dir() +
                                           '/source',
                                           'doc/sphinx/html' ],
                               install : true,
                               install_dir : man_path )

#######################################################################
