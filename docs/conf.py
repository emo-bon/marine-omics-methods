# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
# import os
import sys
from pathlib import Path

# for conversion from markdown to html
# sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, str(Path("..").resolve()))

project = "marine-omics"
copyright = "2025, David Palecek"
author = "David Palecek"
release = "0.1.2"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # 'autoapi.extension',
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]
# autoapi_dirs = ['../momics']

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
