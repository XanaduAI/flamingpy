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

There is a vast literature available to understand the theoretical concepts behind FlamingPy. For a self-contained description, see Xanadu's `blueprint <https://quantum-journal.org/papers/q-2021-02-04-392/>`__ for a fault-tolerant photonic quantum computer.

To see a sample of what FlamingPy can do, let us first import a few important objects:

.. code-block:: python3

   from flamingpy.codes import SurfaceCode
   from flamingpy.cv.ops import CVLayer
   from flamingpy.decoders import correct

Next, let us instantiate an RHG lattice -- the measurement-based version of the surface code:

.. code-block:: python3

   RHG = SurfaceCode(3)


The integer denotes the code distance. By default, the boundaries are set to "open". Next, let us associate the nodes in the RHG lattice with CV states:

.. code-block:: python3

   CVRHG = CVLayer(RHG.graph, p_swap=0.5)


Now, half the lattice (on average) will be labelled a GKP state, and the other half a p-squeezed state. Next, we can apply a noise model to the states:

.. code-block:: python3
   
   grn_model = {"noise": "grn", "delta": 0.1}
   CVRHG.apply_noise(grn_model)


This results in Gaussian random noise model with a squeezing parameter of 0.1 to the GKP states in the lattice. We can now conduct a homodyne measurement on the lattice to measure the syndrome:

.. code-block:: python3

   CVRHG.measure_hom("p", RHG.all_syndrome_inds)


At this point, we are ready to perform error correction on the lattice. First, we can specify some options for the decoder:

.. code-block:: python3

   decoder = {"inner": "basic", "outer": "MWPM"}


This corresponds to a basic GKP binning function for the inner decoder, and minimum-weight perfect matching (MWPM) for the outer decoder. Lastly, we can detect and correct for errors, and print a message identifying success or failure:

.. code-block:: python3

   c = correct(code=RHG, decoder=decoder)
   outcome = "succeeded." * bool(c) + "failed." * (1 - bool(c))
   message = "Error correction {}".format(outcome)
   print(message)


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
   source/fp.cpp.lemonpy
