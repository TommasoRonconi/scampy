.. image:: https://travis-ci.com/TommasoRonconi/scampy.svg?token=snqFgNBsFa2kbkxnNmAN&branch=master
    :target: https://travis-ci.com/TommasoRonconi/scampy

.. image:: https://img.shields.io/badge/ascl-2002.006-blue.svg?colorB=262255
    :target: http://ascl.net/2002.006

.. image:: https://readthedocs.org/projects/scampy/badge/?version=latest
    :target: https://scampy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

A Python interface for Sub-halo Clustering and Abundance Matching
-----------------------------------------------------------------

ScamPy is a highly-optimized and flexible library for "painting" an observed population of cosmological objects on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations.
The method used combines the classical Halo Occupation Distribution (HOD) with the sub-halo abundance matching (SHAM), the sinergy of the two processes is dubbed Sub-halo clustering and
abundance matching (SCAM).
The procedure itself is quite easy since it only requires to apply the two methods in sequence:

1. by applying the HOD scheme, the host sub-haloes are selected;
2. the SHAM algorithm associates to each sub-halo an observable property of choice.

What can be achieved
^^^^^^^^^^^^^^^^^^^^

Here is an animation obtained by running ScamPy on the halo/subhalo catalogues of 42 different snapshots, from redshift z=8 to redshift z=0, of the same :math:`64 Mpc/h` DM-only N-body simulation.
The simulation has been obtained with the non-public code GADGET-3, following the evolution of :math:`512^3` DM particles.
For each different redshift we have fixed the parameters values for the HOD and matched the UV-luminosity function of star-forming galaxies.

.. image:: https://raw.githubusercontent.com/TommasoRonconi/scampy/master/plots/evolving_slice.gif
   :width: 100%
   :alt: If the image is not directly shown in the text, it can be found in the subdirectory `plots/evolving_slice.gif`

The background color-code shows the underlying DM-density field computed by smoothing the contribution of DM-particles in a 10 Mpc/h thick slice of the simulation.
The markers locate the positions of the mock galaxies generated with ScamPy. Circles mark the position of the central galaxies while crosses mark the position of satellite galaxies.
The marker color represents lower to higher luminosity going from brighter to darker.

   
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

You can decide to install it either in
- :code:`developer-mode`, with shared libraries for the C/C++ sectors and headers organized in the
  POSIX directory structure (libraries in :code:`lib`, headers in :code:`include`, python package in :code:`lib/pythonX.Y/site-packages`)
- :code:`package-mode`, with the C++ sector compiled into static libraries within an internal sub-module of the
  package and C-wrapping compiled dynamically along with the former. This is what you would obtain by :code:`pip`-installing from the root
  directory of the project.

.. references:
   
.. _`Meson Build System`: https://mesonbuild.com/

:code:`developer-mode` Install
''''''''''''''''''''''''''''''

From the root directory of this repository run

.. code-block:: bash
		
   meson build_dir --prefix /path/to/install_directory
   meson install -C build_dir

If no :code:`--prefix` is specified the library will be installed in the system default prefix directory (usually :code:`/usr/local`).

:code:`package-mode` Install
''''''''''''''''''''''''''''

From the root directory of this repository either run

.. code-block:: bash
		
   meson build_dir --prefix /path/to/install_directory -Dfull-build=false
   meson install -C build_dir

or run

.. code-block:: bash
   pip install .

In the latter case the standard path for the python installation directory will be used.

Meson options
'''''''''''''

- :code:`full-build`: *boolean*, enables/disables the full build installation.
- :code:`enable-doc`: *boolean*, enables/disables building of the documentation. If enabled, docs will appear in the :code:`$PREFIX/share/man` directory
- :code:`enable-test`: *boolean*, enables/disables testing (to run tests after having compiled the project run :code:`meson test -C build_dir` from the root directory of this repository)

Pre-requisites
''''''''''''''

**For building:**

- :code:`meson<0.57` build system tool
- :code:`ninja`

can both be installed either via :code:`conda install` or with :code:`pip install`

.. warning::
   The current latest version of :code:`meson` (i.e. :code:`0.57.2`) does not always support compiling heritage fortran programs
   (typically an error of type :code:`UnicodeDecodeError` is raised). 
   If the external library FFTLog (see below) is not already installed in your system (and visible to the linker),
   the installation process will try to download and compile it with :code:`ninja`.
   If your :code:`meson` version is superior to :code:`0.56.2` this will cause a failure in the installation process.
   The quickest fix is to downgrade your build system tool to :code:`meson<=0.56.2`.

**Dependencies of the library:**

- GNU Scientific Library version 2 or greater (`GSL link`_); 
- FFTLog (`FFTLog link`_).

While GSL has to be already installed in the system, if FFTLog is not present Meson will authomatically download it along with a patch_ we have developed, both will be installed in the :code:`subprojects` directory of the repository.

**Dependencies for building the documentation locally:**

- Doxygen
- Sphynx (with breathe, autodoc and rtd_theme extensions)

.. note::
  
   A YAML file containing the specs for building a conda environment with all the dependencies needed to build the docs is available at `doc_environment.yml <https://raw.githubusercontent.com/TommasoRonconi/documentation/master/useful/doc_environment.yml>`_

**Dependencies for enabling testing:**

- Google Test (is authomatically installed by Meson)
  
.. references:

.. _patch: https://github.com/TommasoRonconi/fftlog_patch
.. _`GSL link`: https://www.gnu.org/software/gsl/
.. _`FFTLog link`: https://jila.colorado.edu/~ajsh/FFTLog/index.html

