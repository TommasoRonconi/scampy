# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import glob
import os
import subprocess
import sys
from datetime import datetime
# sys.path.insert(0, os.path.abspath('..'))
import scampy
import scampy.io
import scampy.halo
import scampy.utilities
import scampy.measure
import scampy.plot

# -- Configuration for ReadTheDocs setup -------------------------------------

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

# -- Convert tutorials from Jupyter notebooks (only if examples exist) --------

# for fn in glob.glob("../examples/*.ipynb"):
#     name = os.path.splitext(os.path.split(fn)[1])[0]
#     outfn = os.path.join("tutorials", name + ".rst")
#     print("Building {0}...".format(name))
#     subprocess.check_call(
#         "jupyter nbconvert --to rst " + fn + " --output-dir tutorials",
#         shell=True,
#     )
#     subprocess.check_call("python fix_internal_links.py " + outfn, shell=True)

# -- Project information -----------------------------------------------------

project = 'SCAMPy'
copyright = f'{datetime.now().year}, Tommaso Ronconi'
author = 'Tommaso Ronconi'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    # 'sphinx.ext.graphviz'
]

# Autosummary configuration commands:
autodoc_member_order = 'bysource'
autosummary_member_order = 'bysource'

# GraphViz configuration commands:
# graphviz_output_format = 'svg'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = ['.rst', '.md']

master_doc = 'index'

# -- Autodoc configuration ---------------------------------------------------

autodoc_mock_imports = [ 'numpy', 'scipy' ]
autodoc_default_options = {
    'members': True
    }

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

# Define position of the logo
html_logo = '../images/scampy_logo.png'
