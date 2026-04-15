# ScamPy

A Python interface for Sub-halo Clustering and Abundance Matching.

ScamPy is a highly-optimized and flexible library for "painting" an observed population of cosmological objects on top of the DM-halo/subhalo hierarchy obtained from DM-only simulations. The method combines the classical Halo Occupation Distribution (HOD) with sub-halo abundance matching (SHAM) — a synergy dubbed Sub-halo Clustering and Abundance Matching (SCAM).

## Method

1. By applying the HOD scheme, host sub-haloes are selected.
2. The SHAM algorithm associates an observable property of choice to each sub-halo.

## Installation

```bash
pip install pybind11
python setup.py install
```

Or using the provided conda environment:

```bash
conda env create -f scampy_environment.yml
conda activate scampy
python setup.py install
```

## Requirements

- numpy
- scipy
- matplotlib
- h5py
- pybind11 (build-time)
