Documentation
=============

The FlamingPy documentation is built using sphinx.

To build the documentation locally, the following additional packages are required:

* `Sphinx <http://sphinx-doc.org/>`_ >=1.5
* `graphviz <http://graphviz.org/>`_ >=2.38
* `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ >=0.3.6

If using Ubuntu, they can be installed via a combination of ``apt`` and ``pip``:
::

    $ sudo apt install graphviz
    $ sudo apt install python3-sphinx
    $ pip3 install sphinxcontrib-bibtex --user

To build the HTML documentation, go to the `doc` directory and run the command
::

  $ make html

The documentation can then be found in the ``doc/_build/html/`` directory.
