Getting Started and Basic Usage
===============================

There is a vast literature available to understand the theoretical concepts behind FlamingPy. For a self-contained description, see Xanadu's |blueprint| for a fault-tolerant photonic quantum computer.

To see a sample of what FlamingPy can do, let us first import a few important objects:

.. code-block:: python3

   from flamingpy.codes import SurfaceCode
   from flamingpy.decoders import correct
   from flamingpy.noise import CVLayer

Next, let us instantiate an RHG lattice -- the measurement-based version of the surface code:

.. code-block:: python3

   RHG = SurfaceCode(3)


The integer denotes the code distance. By default, the boundaries are set to "open". Next, let us associate the nodes in the RHG lattice with CV states:

.. code-block:: python3

   CVRHG = CVLayer(RHG, delta=0.1, p_swap=0.5)
   CVRHG.apply_noise()


This had the effect of labelling half the lattice (on average) with GKP states and the other half with p-squeezed states. Then, a Gaussian random noise model was applied with a squeezing parameter of 0.1 to the states in the lattice. Finally, a syndrome measurement (sequence of homodyne measurements) was conducted on the lattice, with the outcomes translated to bit values.

At this point, we are ready to perform error correction on the code and print a message identifying success or failure:

.. code-block:: python3
   
   c = correct(RHG)
   outcome = "succeeded." * bool(c) + "failed." * (1 - bool(c))
   message = "Error correction {}".format(outcome)
   print(message)


.. |blueprint| raw:: html

   <a href="https://quantum-journal.org/papers/q-2021-02-04-392/" target="_blank">blueprint</a>