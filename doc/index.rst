.. FlamingPy documentation master file, created by
   sphinx-quickstart on Tue Feb  1 17:32:06 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

FlamingPy
=========


.. image:: _static/RHG_matching.svg
    :align: right
    :width: 250px
    :target: javascript:void(0);

FlamingPy is a cross-platform Python library with several backends for efficient simulations of error correction in fault-tolerant quantum computers

Features
--------

* Simulates error correction on combinations of CV and DV codes to obtain estimations of fault-tolerant thresholds. 
* Supports encoding qubits into GKP states (more precisely, combinations of GKP and squeezed states). 
* Is conveniently modularized, allowing the user to insert custom noise models, codes, decoders, backends and other features. 
* Provides a host of visualization tools for ease of verifying correctness.

Getting started and basic usage
-------------------------------

Coming soon ...

Tutorials
---------

Coming soon ...


.. toctree::
   :maxdepth: 1
   :caption: Using FlamingPy
   :hidden:

   introduction/introduction

.. toctree::
   :maxdepth: 1
   :caption: Development
   :hidden:

   development/development_guide

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
   source/fp.lemonpy
