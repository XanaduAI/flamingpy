# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import subprocess
from unittest import mock

sys.path.insert(0, os.path.abspath("_ext"))
sys.path.insert(0, os.path.abspath("/opt/miniconda/envs/pie/lib/python3.9/site-packages"))
sys.path.append(os.path.abspath(".."))
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(".")), "doc"))


# -- Project information -----------------------------------------------------

project = "FlamingPy"
copyright = "2022, Xanadu Inc."
author = "Xanadu Inc."

# The full version, including alpha/beta/rc tags.
from flamingpy import __version__ as release
# The short X.Y version.
version = re.match(r"^(\d+\.\d+)", release).expand(r"\1")


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.imgmath",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "edit_on_github",
]

MOCK_MODULES = [
    "lemonpy",
]
for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock.Mock()

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "setup.py"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
html_theme = "xanadu_theme"
html_theme_path = ["."]

# Register the theme as an extension to generate a sitemap.xml
# extensions.append("guzzle_sphinx_theme")

# xanadu theme options (see theme.conf for more information)
html_theme_options = {
    # Set the path to a special layout to include for the homepage
    # "homepage": "special_index.html",
    # Set the name of the project to appear in the left sidebar.
    "project_nav_name": "flamingpy",
    "touch_icon": "_static/logo_new.png",
    # colors
    "navigation_button": "#b13a59",
    "navigation_button_hover": "#712b3d",
    "toc_caption": "#b13a59",
    "toc_hover": "#b13a59",
    "table_header_bg": "#ffdce5",
    "table_header_border": "#b13a59",
    "download_button": "#b13a59",
}

edit_on_github_project = "XanaduAI/flamingpy"
edit_on_github_branch = "main/doc"


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
# html_sidebars = {}
html_sidebars = {
    "**": [
        "logo-text.html",
        "searchbox.html",
        "globaltoc.html",
    ]
}

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "FlamingPydoc"


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',
    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',
    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, "FlamingPy.tex", "FlamingPy Documentation", "Xanadu Inc.", "manual"),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [(master_doc, "flamingpy", "FlamingPy Documentation", [author], 1)]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "FlamingPy",
        "FlamingPy Documentation",
        author,
        "FlamingPy",
        "One line description of project.",
        "Miscellaneous",
    ),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ["search.html"]


# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {"https://docs.python.org/": None}

# the order in which autodoc lists the documented members
autodoc_member_order = "bysource"

# inheritance_diagram graphviz attributes
inheritance_node_attrs = dict(color="lightskyblue1", style="filled")

from custom_directives import CustomGalleryItemDirective, DetailsDirective


def setup(app):
    app.add_directive("customgalleryitem", CustomGalleryItemDirective)
    app.add_directive("details", DetailsDirective)
    app.add_css_file("xanadu_gallery.css")
