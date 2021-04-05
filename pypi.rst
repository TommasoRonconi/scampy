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

Here is an animation obtained by running ScamPy on the halo/subhalo catalogues of 42 different snapshots, from redshift z=8 to redshift z=0, of the same `64 Mpc/h` DM-only N-body simulation.
The simulation has been obtained with the non-public code GADGET-3, following the evolution of `512^3` DM particles.
For each different redshift we have fixed the parameters values for the HOD and matched the UV-luminosity function of star-forming galaxies.

.. image:: https://raw.githubusercontent.com/TommasoRonconi/scampy/master/plots/evolving_slice.gif
   :width: 100%

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

Dependencies
^^^^^^^^^^^^

If you are installing from Source Distribution (:code:`sdist`) be aware of the dependencies!

The compiled sector of the package depends on the `GNU Scientific Library (GSL) <https://www.gnu.org/software/gsl/>`_ and thus it should be installed and visible to your system to allow for correct runtime linking.

Precompiled binary packages are included in most GNU/Linux distributions.
A compiled version of GSL is available as part of Cygwin on Windows.

On Linux Debian distributions GSL can be installed using :code:`apt`:

.. code-block:: bash

  $ sudo apt install libgsl-dev -y
