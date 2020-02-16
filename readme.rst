.. image:: https://travis-ci.com/TommasoRonconi/scampy.svg?token=snqFgNBsFa2kbkxnNmAN&branch=master
    :target: https://travis-ci.com/TommasoRonconi/scampy

A Python interface for Sub-halo Clustering and Abundance Matching
-----------------------------------------------------------------

ScamPy is a highly-optimized and flexible library for "painting" an observed population of cosmological objects on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations.
The method used combines the classical Halo Occupation Distribution (HOD) with the sub-halo abundance matching (SHAM), the sinergy of the two processes is dubbed Sub-halo clustering and
abundance matching (SCAM).
The procedure itself is quite easy since it only requires to apply the two methods in sequence:

1. by applying the HOD scheme, the host sub-haloes are selected;
2. the SHAM algorithm associates to each sub-halo an observable property of choice.

Basic Usage
^^^^^^^^^^^

If you wanted to populate a DM catalogue with galaxies with given luminosity, you would do something like:

.. code-block:: python

   # read the sub-halo catalogue from a file
   from scampy import catalogue
   cat = catalogue.catalogue()
   cat.read_hierarchy_from_gadget( "/path/to/input_directory/subhalo_tab_snap" )
   volume = 512.**3 # for a box with side-lenght = 512 Mpc/h


   # build an object of type occupation probability with given parameters
   from scampy import occupation_p
   ocp = occupation_p.tinker10_p( Amin    = 1.e+14,
                                  siglogA = 0.5,
				  Asat    = 1.e+15,
				  alpsat  = 1. )

   # populate the catalogue
   galaxies = cat.populate( ocp, extract = True )

   # define a Schechter-like luminosity function
   import numpy as np
   def schechter ( mag ) :
	alpha = -1.07
	norm = 1.6e-2
	mstar = -19.7 + 5. * np.log10( 5. )
	lum = - 0.4 * ( mag - mstar )
	return 0.4 * np.log( 10 ) * norm * 10**( - 0.07 * lum ) * np.exp( - 10**lum )

   # call the sub-halo abundance matching routine:
   from scampy import abundance_matching
   galaxies = abundance_matching.abundance_matching( galaxies, schechter, factM = 1. / volume )

The :code:`galaxies` array contains the output mock-galaxies.

Installation guide
^^^^^^^^^^^^^^^^^^

Installation of ScamPy is dealt by the  `Meson Build System`_.
Each module of the API is built by a specific :code:`meson.build` script.

.. references:
   
.. _`Meson Build System`: https://mesonbuild.com/


Pre-requisites
''''''''''''''

**For building:**

- :code:`meson` build system tool
- :code:`ninja`

can both be installed either via :code:`conda install` or with :code:`pip install`

**Dependencies of the library:**

- GNU Scientific Library version 2 or greater (`GSL link`_); 
- FFTLog (`FFTLog link`_).

While GSL has to be already installed in the system, if FFTLog is not present Meson will authomatically download it along with a patch_ we have developed, both will be installed in the :code:`subprojects` directory of the repository.

**Dependencies for building the documentation locally:**

- Doxygen
- Sphynx (with breathe, autodoc and rtd_theme extensions)


**Dependencies for enabling testing:**

- Google Test (is authomatically installed by Meson)
  
.. references:

.. _patch: https://github.com/TommasoRonconi/fftlog_patch
.. _`GSL link`: https://www.gnu.org/software/gsl/
.. _`FFTLog link`: https://jila.colorado.edu/~ajsh/FFTLog/index.html

Install
'''''''

.. code-block:: bash
		
   # from root directory
   meson build_dir -Dprefix=/path/to/install_directory
   cd build_dir

   ninja install

If no :code:`-Dprefix` is specified the library will be installed in :code:`/path/to/scampy/install` directory.

