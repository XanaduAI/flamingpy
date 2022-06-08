Build Documentation
===================

The FlamingPy documentation is built using ``sphinx``. To build the documentation locally, the following packages (some in addition to the packages in the top ``dev_requirements.txt``) are required:

* |docutils| >= 0.15.2
* |m2r2| >= 0.3.2
* |matplotlib| >= 3.3.3
* |networkx| >= 2.5
* |numba| >= 0.53.1
* |NumPy| >= 1.21.0
* |pylint| >= 2.13.9
* |retworkx| >= 0.10.2
* |scipy| >= 1.6.0
* |sphinx| >= 4.3.1
* |sphinxcontrib_bibtex| >= 0.4.2
* |sphinx_autodoc_typehints| >= 1.10.3
* |sphinx-automodapi| >= 0.13
* |sphinx-copybutton| >= 0.4
* |thewalrus| >= 0.19.0
* |toml| >= 0.10.2

All required packages can be installed via:
::

    $ python -m pip install -r doc/dev_requirements.txt

We will build the documentation using ``make`` (if you require to install and understand how Makefiles work, see for example |makefile_guide|). Now to build the HTML documentation, go to the ``doc`` directory and run
::

  $ make docs

Note this command will initially run ``make clean`` in the directory.

The documentation can be found in the :file:`doc/_build/html/` directory.


.. |docutils| raw:: html

   <a href="https://docutils.sourceforge.io/" target="_blank">docutils</a>

.. |m2r2| raw:: html

   <a href="https://pypi.org/project/m2r2/" target="_blank">m2r2</a>

.. |matplotlib| raw:: html

   <a href="https://matplotlib.org/" target="_blank">matplotlib</a>

.. |networkx| raw:: html

   <a href="https://networkx.org/" target="_blank">networkx</a>

.. |numba| raw:: html

   <a href="https://numba.pydata.org/" target="_blank">numba</a>

.. |NumPy| raw:: html

   <a href="http://numpy.org/" target="_blank">NumPy</a>

.. |pylint| raw:: html

   <a href="https://pypi.org/project/pylint/" target="_blank">pylint</a>

.. |retworkx| raw:: html

   <a href="https://qiskit.org/documentation/retworkx/" target="_blank">retworkx</a>

.. |scipy| raw:: html

   <a href="https://scipy.org/" target="_blank">scipy</a>

.. |sphinx| raw:: html

   <a href="https://www.sphinx-doc.org/en/master/index.html" target="_blank">sphinx</a>

.. |sphinxcontrib_bibtex| raw:: html

   <a href="https://sphinxcontrib-bibtex.readthedocs.io/en/latest/" target="_blank">sphinxcontrib-bibtex/a>

.. |sphinx_autodoc_typehints| raw:: html

   <a href="https://pypi.org/project/sphinx-autodoc-typehints/" target="_blank">sphinx_autodoc_typehints/a>

.. |sphinx_automodapi| raw:: html

   <a href="https://sphinx-automodapi.readthedocs.io/en/latest/" target="_blank">sphinx-automodapi/a>

.. |sphinx_copybutton| raw:: html

   <a href="https://sphinx-copybutton.readthedocs.io/en/latest/" target="_blank">sphinx-copybutton/a>

.. |thewalrus| raw:: html

   <a href="https://the-walrus.readthedocs.io/en/latest/" target="_blank">thewalrus/a>

.. |toml| raw:: html

   <a href="https://pypi.org/project/toml/" target="_blank">toml/a>

.. |makefile_guide| raw:: html

   <a href="https://pakstech.com/blog/make-windows/#:~:text=make%20%3A%20The%20term%20'make',choose%20Path%20and%20click%20Edit." target="_blank">this guide/a>