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
sys.path.insert(0, os.path.abspath('..'))

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
copyright = '2020, Tommaso Ronconi'
author = 'Tommaso Ronconi'

# -- General configuration ---------------------------------------------------

extensions = [
    'myst_parser',
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = ['.rst', '.md']

master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']
