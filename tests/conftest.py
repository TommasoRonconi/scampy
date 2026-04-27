"""Shared fixtures for scampy regression tests."""

import numpy as np
import pytest

from scampy.cosmology import model as CosmoModel


@pytest.fixture(scope="session")
def default_cosmo():
    """Standard flat ΛCDM cosmology (Ω_M=0.3, Ω_Λ=0.7, h=0.7)."""
    return CosmoModel()


@pytest.fixture(scope="session")
def eds_cosmo():
    """Einstein–de Sitter cosmology (Ω_M=1, Ω_Λ=0, h=1)."""
    return CosmoModel(Om_M=1.0, Om_L=0.0, Om_K=0.0, hh=1.0)


@pytest.fixture(scope="session")
def flat_cosmo_params():
    return dict(Om_M=0.3, Om_b=0.045, Om_L=0.7, hh=0.7, sigma8=0.8)


@pytest.fixture(scope="session")
def simple_pk(default_cosmo):
    """Harrison–Zel'dovich-like power spectrum P(k) ∝ k for basic halo tests."""
    from scampy.power_spectrum import power_spectrum
    kh = np.geomspace(1e-4, 1e2, 500)
    pk0 = kh  # P(k) ∝ k^1, normalized below
    ps = power_spectrum(kh, pk0, cosmo=default_cosmo)
    return ps
