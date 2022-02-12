Development guide
=================

Dependencies
------------

FlamingPy requires the following libraries be installed:

* `Python <http://python.org/>`_ >= 3.8

as well as the following Python packages:

* `NumPy <http://numpy.org/>`_ >= 1.21.0
* `The Walrus <https://the-walrus.readthedocs.io>`_ == 0.15.0
* `matplotlib <https://matplotlib.org/>`_ == 3.3.3
* `networkx <https://networkx.org/>`_ == 2.5
* `retworkx <https://qiskit.org/documentation/retworkx/>`_ == 0.10.2
* `pandas <https://pandas.pydata.org/>`_ == 1.2.1
* `scipy <https://scipy.org/>`_ == 1.6.0


If you currently do not have Python 3 installed, we recommend
`Anaconda for Python 3 <https://www.anaconda.com/download/>`_, a distributed version
of Python packaged for scientific computation.

Development environment
-----------------------

FlamingPy uses a ``pytest`` suite for testing and ``pytest-cov`` for code coverage. These dependencies among compilation 
and other requirements can be installed via ``pip``:

.. code-block:: bash

    pip install -r dev_requirements.txt

Installation
------------

For development purposes, it is recommended to install the FlamingPy source code
using development mode:

.. code-block:: bash

    git clone https://github.com/XanaduAI/ft-stack
    cd ft-stack
    pip install -e .

The ``-e`` flag ensures that edits to the source code will be reflected when
importing FlamingPy in Python.

Software tests
--------------

The FlamingPy test suite includes `pytest <https://docs.pytest.org/en/latest/>`_
and `pytest-cov <https://pytest-cov.readthedocs.io/en/latest/>`_ for coverage reports.

To ensure that FlamingPy is working correctly after installation, the test suite
can be run by navigating to the source code folder and running

.. code-block:: bash

    python3 -m pytest tests


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
``doc/requirements.txt``. These packages can be installed using:

.. code-block:: bash

    pip install -r doc/requirements.txt

from the `doc` directory to then build the HTML documentation, run

.. code-block:: bash

    make html

The documentation can be found in the :file:`doc/_build/html/` directory.
