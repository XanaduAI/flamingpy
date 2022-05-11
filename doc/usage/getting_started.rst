Getting Started and Basic Usage
===============================

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

   CVRHG = CVLayer(RHG, p_swap=0.5)


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
