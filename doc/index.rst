.. FlamingPy documentation master file, created by
   sphinx-quickstart on Tue Feb  1 17:32:06 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FlamingPy
=========

FlamingPy is a cross-platform Python library with a variety of backends for efficient simulations of error correction in fault-tolerant quantum computers

Features
--------

.. image:: _static/RHG_matching.svg
    :align: right
    :width: 250px
    :target: javascript:void(0);

* Simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. 
* Supports encoding qubits into GKP states (more precisely, combinations of GKP and squeezed states). 
* Is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features. 
* Provides a host of visualization tools for ease of verifying correctness.


.. toctree::
   :maxdepth: 1
   :caption: Home
   :hidden:

   self
   easy_installation

.. toctree::
   :maxdepth: 1
   :caption: Using FlamingPy
   :hidden:

   usage/getting_started
   usage/tutorials

.. toctree::
   :maxdepth: 1
   :caption: Development
   :hidden:

   development/guide_for_devs
   build_docs
   development/contribution
   
.. toctree::
   :maxdepth: 1
   :caption: Getting Help
   :hidden:

   help/support
   help/frequently_encountered_errors

.. toctree::
   :maxdepth: 1
   :caption: API
   :hidden:

   source/fp
   source/fp.codes
   source/fp.cv
   source/fp.decoders
   source/fp.simulations
   source/fp.utils
   source/fp.cpp.lemonpy
