.. image:: https://travis-ci.com/TommasoRonconi/scampy.svg?token=snqFgNBsFa2kbkxnNmAN&branch=master
    :target: https://travis-ci.com/TommasoRonconi/scampy

A Python API for Sub-halo Clustering and Abundance Matching
-----------------------------------------------------------

ScamPy is a highly-optimized and flexible library for "painting" an observed population of cosmological objects on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations.
The method used combines the classical Halo Occupation Distribution (HOD) with the sub-halo abundance matching (SHAM), the sinergy of the two processes is dubbed Sub-halo clustering and
abundance matching (SCAM).
The procedure itself is quite easy since it only requires to apply the two methods in sequence:

1. applying the HOD scheme the host sub-haloes are selected;
2. the SHAM algorithm associate to each sub-halo an observable property of choice.

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
