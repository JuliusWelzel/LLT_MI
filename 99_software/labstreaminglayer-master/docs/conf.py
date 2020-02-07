# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------

# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Labstreaminglayer'
copyright = '2019, Christian Kothe, David Medine, Chadwick Boulay, Matthew Grivich, Tristan Stenner'
author = 'Christian Kothe, David Medine, Chadwick Boulay, Matthew Grivich, Tristan Stenner'

# The full version, including alpha/beta/rc tags
release = '1.13'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
#    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

html_theme_options = {
#    'logo': 'logo.png',
    'github_user': 'sccn',
    'github_repo': 'labstreaminglayer',
    'github_button': 'true',
    'extra_nav_links': {
        'C++ API': 'https://labstreaminglayer.readthedocs.io/projects/liblsl/',
        }
    }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

master_doc = 'index'  # for Sphinx < 2.0

# intersphinx
intersphinx_mapping = {
    'liblsl': ('https://labstreaminglayer.readthedocs.io/projects/liblsl', None),
    }

