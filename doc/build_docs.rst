Build Documentation
===================

The FlamingPy documentation is built using ``sphinx``. To build the documentation locally, the following packages (some in addition to the packages in the top ``dev_requirements.txt``) are required:

* `docutils <https://docutils.sourceforge.io/>`_ >= 0.15.2
* `m2r2 <https://pypi.org/project/m2r2/>`_ >= 0.3.2
* `matplotlib <https://matplotlib.org/>`_ >= 3.3.3
* `networkx <https://networkx.org/>`_ >= 2.5
* `numba <https://numba.pydata.org/>`_ >= 0.53.1
* `NumPy <http://numpy.org/>`_ >= 1.21.0
* `pylint <https://pypi.org/project/pylint/>`_ >= 2.13.9
* `retworkx <https://qiskit.org/documentation/retworkx/>`_ >= 0.10.2
* `scipy <https://scipy.org/>`_ >= 1.6.0
* `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ >= 4.3.1
* `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ >= 0.4.2
* `sphinx-autodoc-typehints <https://pypi.org/project/sphinx-autodoc-typehints/>`_ >= 1.10.3
* `sphinx-automodapi <https://sphinx-automodapi.readthedocs.io/en/latest/>`_ >= 0.13
* `sphinx-copybutton <https://sphinx-copybutton.readthedocs.io/en/latest/>`_ >= 0.4
* `thewalrus <https://the-walrus.readthedocs.io/en/latest/>`_ >= 0.19.0
* `toml <https://pypi.org/project/toml/>`_ >= 0.10.2

All required packages can be installed via:
::

    $ python -m pip install -r doc/dev_requirements.txt

We will build the documentation using ``make`` (if you require to install and understand how Makefiles work, see for example `this guide <https://pakstech.com/blog/make-windows/#:~:text=make%20%3A%20The%20term%20'make',choose%20Path%20and%20click%20Edit.>`_). Now to build the HTML documentation, go to the ``doc`` directory and run
::

  $ make html

Note this command will initially run ``make clean`` in the directory.

The documentation can be found in the :file:`doc/_build/html/` directory.
