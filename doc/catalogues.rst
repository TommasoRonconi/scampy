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

   * - Release
     - Paper
     - Catalogue
     - Populations
     - Total sources
     - :math:`z` range
   * - :ref:`RadioMocks <catalogues-radiomocks>`
     - Ronconi et al. (2026)
     - Shallow
     - CoG (AGN + SFG), HIG, cross
     - ~269M
     - CoG: :math:`0 < z \leq 5`; HIG: :math:`0 < z \leq 0.5`
   * - :ref:`RadioMocks <catalogues-radiomocks>`
     - Ronconi et al. (2026)
     - Deep
     - CoG (AGN + SFG), HIG, cross
     - ~1.1B
     - CoG: :math:`0 < z \leq 5`; HIG: :math:`0 < z \leq 0.5`

----

.. _catalogues-radiomocks:

RadioMocks — Ronconi et al. (2026)
------------------------------------

These catalogues are described in Sec. 3.2 of Ronconi et al. (2026, MNRAS, in
press).  They are built on a DEMNUni dark-matter light-cone simulation (Carbone
et al. 2016, Parimbelli et al. 2022) using the full SCAMPy pipeline: a
probabilistic halo-occupation distribution (HOD) selects host sub-haloes, and
sub-halo abundance matching (SHAM) assigns observed radio properties to each
host.

Two catalogues are provided, differing in flux-limit depth:

* the **shallow catalogue** is fully constrained by existing observational data
  (VLA-COSMOS and ALFALFA) and is complete above the flux limits of current
  wide-area surveys;
* the **deep catalogue** extends the calibrated model to the expected sensitivity
  limits of future SKA surveys.

Each catalogue is distributed as a collection of **41 HDF5 files**, one per
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
     - 41 slices, :math:`0 < z \leq 5`
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
     - `[Google Drive — shallow] <PLACEHOLDER_SHALLOW_URL>`_

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
     - 41 slices, :math:`0 < z \leq 5`
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
     - `[Google Drive — deep] <PLACEHOLDER_DEEP_URL>`_

.. _catalogues-radiomocks-format:

Data format
^^^^^^^^^^^

Each redshift-slice file is an HDF5 file.  The table below lists the available
fields.  All positions are in the simulation frame; sky coordinates are derived
from the light-cone geometry described in the paper.

**Continuum catalogue (CoG) fields**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field name
     - Unit
     - Description
   * - ``RA``
     - deg
     - *[placeholder — fill in]*
   * - ``Dec``
     - deg
     - *[placeholder — fill in]*
   * - ``redshift``
     - —
     - *[placeholder — fill in]*
   * - ``S_1p4GHz``
     - Jy
     - *[placeholder — fill in]*
   * - ``population``
     - —
     - *[placeholder — fill in (AGN/SFG flag)]*
   * - ``M_halo``
     - :math:`M_\odot\,h^{-1}`
     - *[placeholder — fill in]*
   * - ``M_subhalo``
     - :math:`M_\odot\,h^{-1}`
     - *[placeholder — fill in]*

**HI catalogue (HIG) fields**

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Field name
     - Unit
     - Description
   * - ``RA``
     - deg
     - *[placeholder — fill in]*
   * - ``Dec``
     - deg
     - *[placeholder — fill in]*
   * - ``redshift``
     - —
     - *[placeholder — fill in]*
   * - ``S_21``
     - Jy·Hz
     - *[placeholder — fill in]*
   * - ``M_HI``
     - :math:`M_\odot\,h^{-1}`
     - *[placeholder — fill in]*
   * - ``M_halo``
     - :math:`M_\odot\,h^{-1}`
     - *[placeholder — fill in]*
   * - ``M_subhalo``
     - :math:`M_\odot\,h^{-1}`
     - *[placeholder — fill in]*
