# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'IHMValidation'
copyright = '2025, Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat'
author = 'Arthur Zalevsky, Sai Ganesan, Benjamin M. Webb, Brinda Vallat'
release = '3.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon', # If using Google/NumPy style docstrings
    'sphinxarg.ext', # Argparse
    'myst_parser',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

import os
import sys
sys.path.insert(0, os.path.abspath('/opt')) # Adjust path as needed
sys.path.insert(0, os.path.abspath('/opt/IHMValidation')) # Adjust path as needed
sys.path.insert(0, os.path.abspath('/home/arthur/work/validation/IHMValidation_3.0_dev')) # Adjust path as needed
sys.path.insert(0, os.path.abspath('/home/arthur/work/validation/IHMValidation_3.0_dev/ihm_validation')) # Adjust path as needed
sys.path.insert(0, os.path.abspath('/home/arthur/work/validation/')) # Adjust path as needed

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

html_context = {
  'display_github': True,
  'github_user': 'salilab',
  'github_repo': 'IHMValidation',
  'github_version': 'branch/dev_3.0/'
}

html_theme_options = {
  'collapse_navigation': False,
}
