# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# -- Project information -----------------------------------------------------
import os
import sys
import pathlib

confpath = pathlib.Path(__file__).parent.resolve()

dt_root = confpath.parent.parent.resolve()

dt_src_root = os.path.join(dt_root,'./bagelfit')

sys.path.insert(0, os.path.abspath(os.path.join(dt_src_root,'src/')))
sys.path.insert(0, os.path.abspath(os.path.join(dt_src_root,'examples/')))
sys.path.insert(0, os.path.abspath(dt_src_root))
sys.path.insert(0, os.path.abspath(dt_root))

project = 'BagelFit'
copyright = '2025, Neelesh Soni'
author = 'Neelesh Soni'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

#extensions = ['sphinx.ext.napoleon', 'sphinx.ext.napoleon',]
extensions = [
'sphinx.ext.autodoc', 
'sphinx.ext.napoleon',
'sphinx.ext.autosummary',
#'sphinx_copybutton',
#'sphinx_toggleprompt',
#'sphinx_pyreverse'
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['../_static']
#html_static_path = []