Guide for Developers
====================

Dependencies
------------

FlamingPy requires the following Python version to be installed:

* |Python| >= 3.8

as well as the following Python packages for development purposes:

* |black| >= 19.3b0
* |cmake| >= 3.14
* |codecov| >= 2.1.12
* |cython| >= 0.29.28
* |docformatter| >= 1.5
* |matplotlib| >= 3.3.3
* |mpi4py| >= 3.1.3 (optional, only for Linux users)
* |networkx| >= 2.5
* |NumPy| >= 1.21
* |pylint| ==2.13.5
* |pytest| >= 6.2
* |pytest-cov| >= 3.0
* |pytest-logger| >= 0.5.1
* |pytest-mock| >= 3.6.1
* |pytest-xdist| >= 1.29.0
* |retworkx| >= 0.10.2
* |setuptools| >= 62.3.2
* |scipy| >= 1.6
* |thewalrus| >= 0.19.0

If you currently do not have Python 3 installed, we recommend
|Anaconda|, a distributed version
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

Using MPI
---------
FamingPy's frontend simulator script, ``simulations.py``, now supports simple and highly-scalable MPI jobs through ``mpi4py``
libraries in a non-intrusive manner. The users who do not have or want MPI, can run ``simulations.py`` single-threaded as
per usual without facing any errors. MPI users can speed up Monte Carlo samplings in EC steps virtually up to as many
processors they can throw at it. The script support jobs both on local machines and large-scale clusters.

To setup FlamingPy's MPI dependencies on a linux system run

.. code-block:: bash

    sudo apt install libopenmpi-dev
    python -m pip install mpi4py>=3.1.3

and continue with the steps described in the following section.

Then, MPI users on their local machines can simply run the following for a 4-processor job:

.. code-block:: bash

    mpirun -np 4 python flamingpy/simulations.py

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
``conda install cmake`` or other means before re-attempting the above. Furthermore,
you may wish to try ``conda install git``. 

Software tests
--------------

The FlamingPy test suite includes |pytest|
and |pytest-cov| for coverage reports.

To ensure that FlamingPy is working correctly after installation, the test suite
can be run by navigating to the source code folder and running

.. code-block:: bash

    python -m pytest tests


Test coverage
^^^^^^^^^^^^^

Test coverage can be checked by running

.. code-block:: bash

    python -m pytest tests --cov=ft_stack --cov-report=xml --cov-report=term-missing -p no:warnings

The output of the above command will show the coverage percentage of each
file, as well as the line numbers of any lines missing test coverage.

To obtain coverage, the ``pytest-cov`` plugin is needed.

Documentation
-------------

See :doc:`build_docs` for the details on how to build the HTML documentation.


.. |Python| raw:: html

   <a href="http://python.org/" target="_blank">Python</a>

.. |black| raw:: html

   <a href="https://pypi.org/project/black/" target="_blank">black</a>

.. |cmake| raw:: html

   <a href="https://pypi.org/project/cmake/" target="_blank">cmake</a>

.. |codecov| raw:: html

   <a href="https://about.codecov.io/language/python/" target="_blank">codecov</a>

.. |cython| raw:: html

   <a href="https://cython.org/" target="_blank">cython</a>

.. |docformatter| raw:: html

   <a href="https://pypi.org/project/docformatter/" target="_blank">docformatter</a>

.. |matplotlib| raw:: html

   <a href="https://matplotlib.org/" target="_blank">matplotlib</a>

.. |mpi4py| raw:: html

   <a href="https://mpi4py.readthedocs.io/en/stable/" target="_blank">mpi4py</a>

.. |networkx| raw:: html

   <a href="https://networkx.org/" target="_blank">networkx</a>

.. |NumPy| raw:: html

   <a href="http://numpy.org/" target="_blank">NumPy</a>

.. |pylint| raw:: html

   <a href="https://pylint.pycqa.org/en/latest/" target="_blank">pytest</a>

.. |pytest| raw:: html

   <a href="https://docs.pytest.org/en/latest/" target="_blank">pytest</a>

.. |pytest-cov| raw:: html

   <a href="https://pypi.org/project/pytest-cov/" target="_blank">pytest-cov</a>

.. |pytest-logger| raw:: html

   <a href="https://pypi.org/project/pytest-logger/" target="_blank">pytest-logger</a>

.. |pytest-mock| raw:: html

   <a href="https://pypi.org/project/pytest-mock/" target="_blank">pytest-mock</a>

.. |pytest-xdist| raw:: html

   <a href="https://pypi.org/project/pytest-xdist/" target="_blank">pytest-xdist</a>

.. |retworkx| raw:: html

   <a href="https://qiskit.org/documentation/retworkx/" target="_blank">retworkx</a>

.. |setuptools| raw:: html

   <a href="https://pypi.org/project/setuptools/" target="_blank">setuptools</a>

.. |scipy| raw:: html

   <a href="https://scipy.org/" target="_blank">scipy</a>

.. |thewalrus| raw:: html

   <a href="https://the-walrus.readthedocs.io/en/latest/" target="_blank">thewalrus</a>

.. |Anaconda| raw:: html

   <a href="https://www.anaconda.com/download/" target="_blank">Anaconda for Python 3</a>
