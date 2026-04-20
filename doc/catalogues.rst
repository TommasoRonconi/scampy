.. _catalogues:

Available Published Catalogues
===============================

This page lists mock galaxy catalogues generated with SCAMPy and made available
for public use.  Catalogues are distributed as collections of HDF5 files, one
per redshift slice, hosted on Google Drive.  Access is granted on request: click
the download link for the catalogue of interest and request access through the
Google Drive interface.

If you use any of these catalogues in your work, please cite the corresponding
paper listed in the table below.

.. note::

   Access requests are reviewed manually.  This allows the authors to keep track
   of the user community and to notify users of updates or corrections.

----

.. _catalogues-summary:

Summary
-------

.. list-table::
   :header-rows: 1
   :widths: 20 20 15 15 15 15

   * - Paper
     - Catalogue
     - Populations
     - Total sources
     - :math:`z` range
   * - `[Ronconi et al. (2026)] https://ui.adsabs.harvard.edu/abs/2026arXiv260325650R/abstract`_
     - Shallow
     - CoG (AGN + SFG), HIG, cross
     - ~269M
     - CoG: :math:`0 < z \leq 5`; HIG: :math:`0 < z \leq 0.5`
   * - `[Ronconi et al. (2026)] https://ui.adsabs.harvard.edu/abs/2026arXiv260325650R/abstract`_
     - Deep
     - CoG (AGN + SFG), HIG, cross
     - ~1.1B
     - CoG: :math:`0 < z \leq 5`; HIG: :math:`0 < z \leq 0.5`

----

.. _catalogues-radiomocks:

RadioMocks — Ronconi et al. (2026)
-----------------------------------

These catalogues are described in Sec. 3.2 of `[Ronconi et al. (2026)] https://ui.adsabs.harvard.edu/abs/2026arXiv260325650R/abstract`_.
They are built on a DEMNUni dark-matter light-cone simulation (Carbone et al. 2016, Parimbelli et al. 2022) using the full SCAMPy pipeline: a
probabilistic halo-occupation distribution (HOD) selects host sub-haloes, and
sub-halo abundance matching (SHAM) assigns observed radio properties to each
host.

Two catalogues are provided, differing in flux-limit depth:

* the **shallow catalogue** is fully constrained by existing observational data
  (VLA-COSMOS and ALFALFA) and is complete above the flux limits of current
  wide-area surveys;
* the **deep catalogue** extends the calibrated model to the expected sensitivity
  limits of future SKA surveys.

Each catalogue is distributed as a collection of **42 HDF5 files**, one per
redshift slice, spanning :math:`0 < z \leq 5` for the continuum population and
:math:`0 < z \leq 0.5` for the HI population.

.. _catalogues-radiomocks-shallow:

Shallow catalogue
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0
   :widths: 35 65

   * - **Paper**
     - Ronconi et al. (2026)
   * - **Simulation**
     - DEMNUni light-cone (Carbone et al. 2016)
   * - **Redshift slices**
     - 42 slices, :math:`0 < z \leq 5`
   * - **Flux limit — continuum**
     - :math:`S^{\lim}_{1.4\,\mathrm{GHz}} \simeq 8 \times 10^{-5}` Jy
   * - **Flux limit — HI line**
     - :math:`S^{\lim}_{21} \simeq 2` Jy·Hz
   * - **Continuum sources (CoG,** :math:`z \leq 5` **)**
     - 217 143 484 (AGN: 30 553 719; SFG: 186 589 765)
   * - **HI sources (HIG,** :math:`z \leq 0.5` **)**
     - 51 509 884
   * - **Cross-catalogue (CoG** :math:`\times` **HIG)**
     - 7 084 155
   * - **Download**
     - `[Google Drive — shallow] https://drive.google.com/drive/folders/17UlI0uvlIjzX0TwgDX2A1A7ECgF9FZI-?usp=sharing`_

.. _catalogues-radiomocks-deep:

Deep catalogue
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 0
   :widths: 35 65

   * - **Paper**
     - Ronconi et al. (2026)
   * - **Simulation**
     - DEMNUni light-cone (Carbone et al. 2016)
   * - **Redshift slices**
     - 42 slices, :math:`0 < z \leq 5`
   * - **Flux limit — continuum**
     - :math:`S^{\lim}_{1.4\,\mathrm{GHz}} \simeq 4 \times 10^{-5}` Jy
   * - **Flux limit — HI line**
     - :math:`S^{\lim}_{21} \simeq 0.3` Jy·Hz
   * - **Continuum sources (CoG,** :math:`z \leq 5` **)**
     - 365 527 186 (AGN: 40 607 134; SFG: 324 920 052)
   * - **HI sources (HIG,** :math:`z \leq 0.5` **)**
     - 719 501 958
   * - **Cross-catalogue (CoG** :math:`\times` **HIG)**
     - 32 531 743
   * - **Download**
     - `[Google Drive — deep] https://drive.google.com/drive/folders/1FouNVty4qduCcBbZRg-YX9oudfGskDr1?usp=sharing`_

.. _catalogues-radiomocks-format:

Data format
^^^^^^^^^^^

Each redshift-slice file is an HDF5 file containing a single group ``galaxies/``
with one dataset per field.  All entries in a slice file represent objects
detected in at least one catalogue (CoG or HIG); the boolean flags ``catCoG``
and ``catHIG`` indicate membership.  Comoving coordinates are in the simulation
frame; sky coordinates are derived from the light-cone geometry described in the
paper.

The files can be read with :func:`scampy.io.hdf5.load_from_hdf5`, which returns
a nested dictionary keyed first by group name, then by field name:

.. code-block:: python

   from scampy.io.hdf5 import load_from_hdf5

   data = load_from_hdf5("Catalogue_FluxLim_z0.01.hdf5")
   galaxies = data["galaxies"]

   # sky coordinates (all objects)
   ra, dec, z = galaxies["RA"], galaxies["Dec"], galaxies["Red"]

   # select the continuum sub-catalogue
   cog_mask = galaxies["catCoG"]
   flux_1400 = galaxies["I1400"][cog_mask]   # µJy

   # select the HI sub-catalogue
   hig_mask = galaxies["catHIG"]
   hi_flux = galaxies["HI flux"][hig_mask]   # mJy·Hz

.. list-table::
   :header-rows: 1
   :widths: 20 12 15 53

   * - Field
     - dtype
     - Unit
     - Description
   * - ``RA``
     - float64
     - deg
     - Right ascension
   * - ``Dec``
     - float64
     - deg
     - Declination
   * - ``Red``
     - float64
     - —
     - Redshift
   * - ``X``
     - float64
     - Mpc/h
     - Comoving x-coordinate
   * - ``Y``
     - float64
     - Mpc/h
     - Comoving y-coordinate
   * - ``Z``
     - float64
     - Mpc/h
     - Comoving z-coordinate
   * - ``Mhalo``
     - float64
     - :math:`M_\odot\,h^{-1}`
     - Host halo mass
   * - ``Msubh``
     - float64
     - :math:`M_\odot\,h^{-1}`
     - Host subhalo mass
   * - ``Parent``
     - int64
     - —
     - Parent halo index
   * - ``IDgadget``
     - int64
     - —
     - GADGET particle ID
   * - ``snapID``
     - int64
     - —
     - Snapshot ID
   * - ``Sect``
     - int64
     - —
     - Sky sector index
   * - ``OptClass``
     - int64
     - —
     - Optical source classification
   * - ``logSFR``
     - float32
     - :math:`\log(M_\odot\,\mathrm{yr}^{-1})`
     - Log star-formation rate
   * - ``catCoG``
     - bool
     - —
     - Member of the continuum catalogue (CoG)
   * - ``catHIG``
     - bool
     - —
     - Member of the HI catalogue (HIG)
   * - ``cenCoG``
     - bool
     - —
     - Central galaxy in CoG
   * - ``cenHIG``
     - bool
     - —
     - Central galaxy in HIG
   * - ``satCoG``
     - bool
     - —
     - Satellite galaxy in CoG
   * - ``satHIG``
     - bool
     - —
     - Satellite galaxy in HIG
   * - ``I1400`` *(CoG)*
     - float32
     - :math:`\mu\mathrm{Jy}`
     - Flux density at 1.4 GHz
   * - ``I150`` *(CoG)*
     - float32
     - :math:`\mu\mathrm{Jy}`
     - Flux density at 150 MHz
   * - ``I3000`` *(CoG)*
     - float32
     - :math:`\mu\mathrm{Jy}`
     - Flux density at 3 GHz
   * - ``L1400`` *(CoG)*
     - float32
     - :math:`\log(\mathrm{W\,Hz}^{-1})`
     - Radio luminosity at 1.4 GHz
   * - ``RadioClass`` *(CoG)*
     - int64
     - —
     - Radio source classification (AGN/SFG)
   * - ``size`` *(CoG)*
     - float32
     - arcsec
     - Angular size
   * - ``trecs_CoG_ID`` *(CoG)*
     - int64
     - —
     - T-RECS continuum source ID
   * - ``HI flux`` *(HIG)*
     - float32
     - :math:`\mathrm{mJy\cdot Hz}`
     - HI line flux
   * - ``HI size`` *(HIG)*
     - float32
     - arcsec
     - HI disc angular size
   * - ``MHI`` *(HIG)*
     - float32
     - :math:`\log(M_\odot)`
     - HI mass
   * - ``MHI_pred`` *(HIG)*
     - float32
     - :math:`\log(M_\odot)`
     - Predicted HI mass (from SHAM)
   * - ``w50`` *(HIG)*
     - float32
     - km/s
     - Line width at 50% peak flux
   * - ``v_max`` *(HIG)*
     - float32
     - km/s
     - Maximum rotation velocity
   * - ``inclination`` *(HIG)*
     - float32
     - deg
     - Disc inclination angle
   * - ``trecs_HIG_ID`` *(HIG)*
     - int64
     - —
     - T-RECS HI source ID
