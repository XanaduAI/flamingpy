Build Documentation
===================

The FlamingPy documentation is built using `sphinx`. To build the documentation locally, the following packages are required:

* `Sphinx <http://sphinx-doc.org/>`_ >=1.5
* `graphviz <http://graphviz.org/>`_ >=2.38
* `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ >=0.3.6

All required packages can be installed via:
::

    $ python -m pip install -r doc/dev_requirements.txt

We will build the documentation using `make` (if you require to install and understand how Makefiles work, see for example `this guide <https://pakstech.com/blog/make-windows/#:~:text=make%20%3A%20The%20term%20'make',choose%20Path%20and%20click%20Edit.>`_). Now to build the HTML documentation, go to the `doc` directory and run the command
::

  $ make html

The documentation can then be found in the ``doc/_build/html/`` directory.
