Sub-halo Clustering and Abundance matching Python API
-----------------------------------------------------

Installation guide
^^^^^^^^^^^^^^^^^^

Pre-requisites
''''''''''''''

- :code:`meson` build system tool
- :code:`ninja`

can both be installed either via :code:`conda install` or with :code:`pip install`

Install
'''''''

.. code-block:: bash
		
   # from root directory
   meson build_dir -Dprefix=/path/to/install_directory
   cd build_dir

   ninja install

If no :code:`-Dprefix` is specified the library will be installed in :code:`/path/to/scampy/install` directory.
