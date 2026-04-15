.. image:: https://img.shields.io/badge/ascl-2002.006-blue.svg?colorB=262255
    :target: http://ascl.net/2002.006

.. image:: https://readthedocs.org/projects/scampy/badge/?version=latest
    :target: https://scampy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

A Python interface for Sub-halo Clustering and Abundance Matching
-----------------------------------------------------------------

SCAMPy is a flexible library for "painting" an observed population of cosmological objects
on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations.
The method combines the classical Halo Occupation Distribution (HOD) with sub-halo
abundance matching (SHAM): a synergy dubbed Sub-halo Clustering and Abundance Matching (SCAM).

The procedure itself is quite easy since it only requires to apply the two methods in sequence:

1. by applying the HOD scheme, the host sub-haloes are selected;
2. the SHAM algorithm associates to each sub-halo an observable property of choice.

What can be achieved
^^^^^^^^^^^^^^^^^^^^

Here is an animation obtained by running SCAMPy on the halo/subhalo catalogues of 42 different
snapshots, from redshift :math:`z=8` to redshift :math:`z=0`, of the same :math:`64\,\mathrm{Mpc}/h`
DM-only N-body simulation.
The simulation has been obtained with the non-public code GADGET-3, following the evolution of
:math:`512^3` DM particles.
For each different redshift we have fixed the parameter values for the HOD and matched the
UV-luminosity function of star-forming galaxies.

.. image:: https://raw.githubusercontent.com/TommasoRonconi/scampy/main/images/evolving_slice.gif
   :width: 100%
   :alt: If the image is not directly shown in the text, it can be found in the subdirectory `images/evolving_slice.gif`

The background colour-code shows the underlying DM-density field computed by smoothing the
contribution of DM-particles in a 10 Mpc/h thick slice of the simulation.
The markers locate the positions of the mock galaxies generated with SCAMPy.
Circles mark the position of central galaxies while crosses mark the position of satellite galaxies.
The marker colour represents lower to higher luminosity going from brighter to darker.

Basic Usage
^^^^^^^^^^^

To populate a DM catalogue with mock galaxies assigned a UV magnitude you would do something like:

.. code-block:: python

   # read the GADGET SUBFIND tables
   from scampy.io.subfind_table import subfind_table
   tab = subfind_table( "/path/to/subhalo_tab_snap" )
   tab.read_all_files()

   # build halo and sub-halo catalogue objects
   import numpy as np
   from scampy.catalogue import haloCat, subhaloCat, catalogue

   haloes = haloCat(
       len( tab.content['GroupMass'] ),
       Mhalo    = tab.content['GroupMass'],
       Rhalo    = tab.content['GroupRadiusM200'],
       X        = tab.content['GroupPos'][:, 0],
       Y        = tab.content['GroupPos'][:, 1],
       Z        = tab.content['GroupPos'][:, 2],
       firstSub = tab.content['GroupFirstSub'],
       numSubs  = tab.content['GroupNsubs'],
   )
   subhaloes = subhaloCat(
       len( tab.content['SubhaloMass'] ),
       Msubh  = tab.content['SubhaloMass'],
       Parent = tab.content['SubhaloGrNr'],
       X      = tab.content['SubhaloPos'][:, 0],
       Y      = tab.content['SubhaloPos'][:, 1],
       Z      = tab.content['SubhaloPos'][:, 2],
   )
   cat = catalogue( haloes, subhaloes, Lbox = 512. )  # Lbox in Mpc/h

   # select host sub-haloes with an HOD
   from scampy.hod import HOD, get_hosts
   hod   = HOD( Mmin=1.e12, siglogM=0.5, M0=1.e12, M1=1.e13, alpha=1. )
   hosts = get_hosts( cat, Pcen=hod.Pcen, Psat=hod.Psat )

   # define a Schechter-like UV luminosity function
   def schechter ( mag ) :
       alpha = -1.07
       norm  = 1.6e-2
       mstar = -19.7 + 5. * np.log10( 5. )
       lum   = -0.4 * ( mag - mstar )
       return 0.4 * np.log( 10 ) * norm * 10**( -0.07 * lum ) * np.exp( -10**lum )

   # assign magnitudes via sub-halo abundance matching
   from scampy.sham import match_distribution
   magnitudes = match_distribution(
       cat.subhaloes['Msubh'][hosts], schechter,
       minF = -25, maxF = -10, factF = 1. / cat.volume,
   )

The :code:`magnitudes` array contains the UV magnitude assigned to each mock galaxy.

Installation
^^^^^^^^^^^^

The recommended installation uses :code:`pip` after cloning the repository:

.. code-block:: bash

   git clone https://github.com/TommasoRonconi/scampy.git
   cd scampy
   pip install pybind11
   python setup.py install

Alternatively, a :code:`conda` environment file is provided:

.. code-block:: bash

   conda env create -f scampy_environment.yml
   conda activate scampy
   python setup.py install

Requirements
^^^^^^^^^^^^

**Runtime:**

* :code:`numpy`
* :code:`scipy`
* :code:`matplotlib`
* :code:`h5py`

**Build-time:**

* :code:`pybind11`

**Documentation (optional):**

* :code:`sphinx>=5.0`
* :code:`sphinx-rtd-theme`
* :code:`myst-parser`
