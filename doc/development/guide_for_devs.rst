Guide for developers
====================

Dependencies
------------

FlamingPy requires the following Python version to be installed:

* `Python <http://python.org/>`_ >= 3.8

as well as the following Python packages:

* `docutils <https://docutils.sourceforge.io/>`_ >= 0.15.2
* `m2r2 <https://pypi.org/project/m2r2/>`_ >= 0.3.2
* `matplotlib <https://matplotlib.org/>`_ >= 3.3.3
* `numba <https://numba.pydata.org/>`_ >= 0.53.1
* `networkx <https://networkx.org/>`_ >= 2.5
* `NumPy <http://numpy.org/>`_ >= 1.21.0
* `pandas <https://pandas.pydata.org/>`_ >= 1.2.1
* `retworkx <https://qiskit.org/documentation/retworkx/>`_ >= 0.10.2
* `scipy <https://scipy.org/>`_ >= 1.6.0
* `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ >= 4.3.1
* `sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/>`_ >= 0.4.2
* `sphinx-autodoc-typehints <https://pypi.org/project/sphinx-autodoc-typehints/>`_ >= 1.10.3
* `sphinx-automodapi <https://sphinx-automodapi.readthedocs.io/en/latest/>`_ >= 0.13
* `sphinx-copybutton <https://sphinx-copybutton.readthedocs.io/en/latest/>`_ >= 0.4
* `toml <https://pypi.org/project/toml/>`_ >= 0.10.2

If you currently do not have Python 3 installed, we recommend
`Anaconda for Python 3 <https://www.anaconda.com/download/>`_, a distributed version
of Python packaged for scientific computation.

Setting up a development environment
------------------------------------

If you are a developer and wish to manipulate and test FlamingPy source code, you need 
to set up a development environment. First, clone the FlamingPy repository. 
Then, create and activate a new virtual environment (if you prefer using an existing 
environment, you may need to uninstall existing FlamingPy builds). If you use **Conda**, 
for example, you may run the following:

.. code-block:: bash

    conda create -n flamingpy python=3.8
    conda activate flamingpy

FlamingPy uses a ``pytest`` suite for testing and ``pytest-cov`` for code coverage. These dependencies among compilation 
and other requirements as stated in the above section can be installed via ``pip``:

.. code-block:: bash

    python -m pip install -r dev_requirements.txt

Installation from Source
------------------------

To install FlamingPy from Source, next change to the directory where FlamingPy was cloned and run:

.. code-block:: bash

    python setup.py develop # only installs Python libraries
    python setup.py build_cython --inplace # [OPTIONAL] compiles Cython-based backends
    python setup.py build_cmake --inplace # [OPTIONAL] compiles CMake-based backends 

The purpose of the commands is as follows:

* The first command installs dependencies for building the project and testing purposes, and can be skipped if already satisfied. 
* The second command (develop) installs FlamingPy Python libraries without compiling the optional backends. 
* The next optional commands compile various FlamingPy backends as required (given you have appropriate compilers pre-installed). 

If you encountered a CMake error, you may need to (re-)install it through 
`conda install cmake` or other means before re-attempting the above. Furthermore, 
you may wish to try `conda install git`. For more detailed instructions and 
recommendations, including how to configure your environments, compilers and 
resolve errors, see our Frequently Encountered Errors page.

Software tests
--------------

The FlamingPy test suite includes `pytest <https://docs.pytest.org/en/latest/>`_
and `pytest-cov <https://pytest-cov.readthedocs.io/en/latest/>`_ for coverage reports.

To ensure that FlamingPy is working correctly after installation, the test suite
can be run by navigating to the source code folder and running

.. code-block:: bash

    python -m pytest tests


Test coverage
^^^^^^^^^^^^^

Test coverage can be checked by running

.. code-block:: bash

    python3 -m pytest tests --cov=ft_stack --cov-report=xml --cov-report=term-missing -p no:warnings

The output of the above command will show the coverage percentage of each
file, as well as the line numbers of any lines missing test coverage.

To obtain coverage, the ``pytest-cov`` plugin is needed.

Documentation
-------------

Additional packages are required to build the documentation, as specified in
``doc/dev_requirements.txt``. These packages can be installed using:

.. code-block:: bash

    pip install -r doc/dev_requirements.txt

from the `doc` directory to then build the HTML documentation, run

.. code-block:: bash

    make clean html

The documentation can be found in the :file:`doc/_build/html/` directory.
