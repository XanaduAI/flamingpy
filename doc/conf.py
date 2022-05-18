# -*- coding: utf-8 -*-
#
# flamingpy configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config
"""
flamingpy configuration file for the Sphinx documentation builder.
"""
import os, sys, re, time
from unittest.mock import MagicMock


# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath("_ext"))
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(".")), "doc"))


class Mock(MagicMock):
    """An auxiliary class to create mocked modules"""

    __name__ = "foo"

    @classmethod
    def __getattr__(cls, name):
        return MagicMock()


# MOCK_MODULES = ["flamingpy.cpp.lemonpy", "flamingpy.cpp.cpp_mc_loop"]
MOCK_MODULES = []
mock = Mock()

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = mock


# -- Project information -----------------------------------------------------

project = "FlamingPy"
copyright = "2022, Xanadu Quantum Technologies"
author = "Xanadu Inc."

# The full version, including alpha/beta/rc tags.
with open("../flamingpy/_version.py") as f:
    release = f.readlines()[-1].split()[-1].strip("\"'")

# The short X.Y version.
version = re.match(r"^(\d+\.\d+)", release).expand(r"\1")


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "3.0"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "edit_on_github",
    "m2r2",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.graphviz",
    "sphinx.ext.imgmath",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx.ext.inheritance_diagram",
    "sphinx_gallery.gen_gallery",
]

intersphinx_mapping = {"https://flamingpy.readthedocs.io/en/latest/": None}

sphinx_gallery_conf = {
    # path to your example scripts
    "examples_dirs": ["tutorials"],
    # path where to save gallery generated examples
    "gallery_dirs": ["tutorials/_out"],
    # execute files that match the following filename pattern,
    # and skip those that don't. If the following option is not provided,
    # all example scripts in the 'examples_dirs' folder will be skiped.
    "filename_pattern": r"/run_",
    "pypandoc": True,
    # first notebook cell in generated Jupyter notebooks
    "first_notebook_cell": (
        "# This cell is added by sphinx-gallery\n"
        "# It can be customized to whatever you like\n"
        "%matplotlib inline"
    ),
    # thumbnail size
    "thumbnail_size": (400, 400),
    'reference_url': {
         # The module you locally document uses None
        'flamingpy': "https://flamingpy.readthedocs.io/en/latest/",
    },
    'backreferences_dir'  : 'backreferences',
    'doc_module'          : ('flamingpy'),
    'junit': '../test-results/sphinx-gallery/junit.xml',
    'download_all_examples': False,
}

automodapi_toctreedirnm = "source/api"
automodsumm_inherited_members = True
autosummary_generate = True
autosummary_imported_members = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
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

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = True

# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/favicon.ico"

# Creates UML diagrams (svg). These are later used in source/fp.rst .   
os.system('pyreverse -o svg -p flamingpy ../flamingpy -d _static --colorized --max-color-depth 1 -k')
time.sleep(0.5)

#changes color of uml diagrams to match the common 'blue' theme
with open('_static/packages_flamingpy.svg', 'r') as file :
  filedata = file.read()
filedata = filedata.replace('aliceblue', '#bde0ff')
filedata = filedata.replace('green', 'black')
with open('_static/packages_flamingpy.svg', 'w') as file:
  file.write(filedata)

with open('_static/classes_flamingpy.svg', 'r') as file :
  filedata = file.read()
filedata = filedata.replace('aliceblue', '#bde0ff')
filedata = filedata.replace('green', 'black')
with open('_static/classes_flamingpy.svg', 'w') as file:
  file.write(filedata)

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "FlamingPydoc"


# -- Options for other HTML outputs ------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "xanadu"

# Xanadu theme options (see theme.conf for more information).
html_theme_options = {
    "navbar_name": "FlamingPy",
    "navbar_logo_colour": "#f57c00",
    "navbar_left_links": [
        {
            "name": "Quantum Error Correction",
            "href": "quantum_error_correction.html",
        },
        {
            "name": "Install",
            "href": "install.html",
        },
        {
            "name": "Documentation",
            "href": "index.html",
            "active": True,
        },
    ],
    "navbar_right_links": [
        {
            "name": "FAQ",
            "href": "faq.html",
            "icon": "fas fa-question",
        },
        {
            "name": "Support",
            "href": "https://github.com/XanaduAI/flamingpy/issues",
            "icon": "fab fa-discourse",
        },
        {
            "name": "GitHub",
            "href": "https://github.com/XanaduAI/flamingpy",
            "icon": "fab fa-github",
        },
    ],
    "prev_next_button_colour": "#f57c00",
    "prev_next_button_hover_colour": "#bb4d00",
    "toc_marker_colour": "#f57c00",
    "table_header_background_colour": "#ffdce5",
    "border_colour": "#f57c00",
    "code_colour": "#ef6c00",
    "text_accent_colour": "#f57c00",
    "extra_copyrights": [
        "TensorFlow, the TensorFlow logo, and any related marks are trademarks of Google Inc."
    ],
}

edit_on_github_project = "XanaduAI/flamingpy"
edit_on_github_branch = "main/doc"

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
html_sidebars = {
    "**": [
        "searchbox.html",
        "globaltoc.html",
    ]
}


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
        "FlamingPy is a cross-platform Python library with a variety of backends for efficient simulations of error correction in fault-tolerant quantum computers.",
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


# -- Options for intersphinx extension ---------------------------------------

# the order in which autodoc lists the documented members
autodoc_member_order = "bysource"

# inheritance_diagram graphviz attributes
inheritance_node_attrs = dict(color="lightskyblue1", fillcolor="lightskyblue1", style="filled")
