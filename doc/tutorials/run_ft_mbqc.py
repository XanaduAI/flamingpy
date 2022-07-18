r""".. _ft_mbqc-tutorial:

Fault-tolerant measurement-based quantum computation
=======================================================

*Author: Joost Bus. Posted: Day Month 2022. Last updated: Day Month 2022.*

"""

##############################################################################
# Quantum error correction and fault tolerance
# --------------------------------------------
#
# To mitigate the physical errors that can (and will) happen during a quantum computation we
# require some kind of error correction. Error correction is a technique of detecting errors and
# reconstructing the logical data without losing any information. It is not exclusive to quantum computing;
# it is also ubiquitous in `"classical" computing <https://www.youtube.com/watch?v=AaZ_RSt0KP8>`_
# and communication. However, it is a stringent requisite in the quantum realm as the systems one
# works with are much more precarious and therefore prone to environmental factors, causing errors.
#
# Due to the peculiarities of quantum physics, we have to be careful when though. First of all, we can
# not simply look inside our quantum computer and see if an error occured. This would collapse the
# wavefunction which caries valuable information. Secondly, we can not make copies of a quantum
# state to create redundancy. This is because of the no-cloning theorem. A whole research field is devoted
# to combat these challenges since Peter Shor's published a seminal paper in 1995 [#ShorQEC1995]_. A
# full coverage of this topic is beyond the scope of this tutorial, but a good place to start is
# `Daniel Gottesman's thesis <https://arxiv.org/abs/quant-ph/9705052>`_ or `this blog post by
# Arthur Pesah <https://arthurpesah.me/blog/2022-01-25-intro-qec-1/>`_ for a more compact
# introduction. Instead, what we will do here is showing how to implement error correction in the
# MBQC framework.
#
# In the measurement-based picture, quantum error correction requires a 3-dimensional cluster state
# [#XanaduBlueprint]_. The error correcting code that you want to implement dictates the structure
# of the cluster state. Let's see how we can implement the famous surface code [#FowlerSurfaceCode]_ [#GoogleQEC2022]_ as
# an example. The cluster state that is associated with this code is known as the the RHG lattice,
# named after its architects Raussendorf, Harrington, and Goyal. We can visualize this cluster
# state with FlamingPy.
#

from flamingpy.codes import SurfaceCode
import matplotlib.pyplot as plt

code_distance = 3
RHG = SurfaceCode(code_distance)
# RHG.draw(backend="plotly") #TODO: uncomment this line after merging #103

##############################################################################
#
# .. raw:: html
#    :file: ../tutorials/ft_mbqc/rhg-graph.html
#

##############################################################################
#
# The actual computation is done by performing single-qubit measurements, as illustrated below. At
# each timestep, we measure all the qubits on one sheet of the lattice. The binary outcomes of these
# measurements determine the measurement bases for future measurements, and the last sheet of the
# lattice encodes the result of the computation which can be read out by yet another measurement!
#
#
# .. figure:: ../tutorials/ft_mbqc/gif_measuring.gif
#    :align: center
#    :width: 75%
#
#    ..
#
#    Performing a computation with measurements using the RHG lattice. [#XanaduBlueprint]_
#


##############################################################################
# Xanadu's approach
# --------------------------------
# Xanadu's path towards a fault-tolerant quantum computer is via a measurement-based scheme with a
# 3-dimensional cluster state using photonics. The main ideas of the architecture are presented in
# [#XanaduBlueprint]_ and the corresponding cluster state is shown in the figure below. One
# interesting aspect of this architecture is the use of a hybrid cluster: a combination of GKP
# states and squeezed states. If you want to dive into more depth, you can have a look at `this
# blog post <https://medium.com/xanaduai/from-a-state-of-light-to-state-of-the-art-the-photonic-path-to-millions-of-qubits-c0e08ca1cb21>`_
# or `this video <https://www.youtube.com/watch?v=SD6TH7GZ1rM>`_.
#
# .. figure:: ../tutorials/ft_mbqc/mbqc_blueprint_full.png
#    :align: center
#    :width: 75%
#
#    ..
#
#
#    Hybrid cluster state proposed in [#XanaduBlueprint]_.
#
#
#

##############################################################################
# References
# ----------
#
#
# .. [#OneWay2001]
#
#     Robert Raussendorf and Hans J. Briegel (2001) *A One-Way Quantum Computer*,
#     `Phys. Rev. Lett. 86, 5188
#     <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.5188>`__.
#
# .. [#XanaduBlueprint]
#
#     J. Eli Bourassa, Rafael N. Alexander, Michael Vasmer et al. (2021) *Blueprint for a Scalable Photonic Fault-Tolerant Quantum Computer*,
#     `Quantum 5, 392
#     <https://quantum-journal.org/papers/q-2021-02-04-392/>`__.
#
# .. [#XanaduPassiveArchitecture]
#
#     Ilan Tzitrin, Takaya Matsuura, Rafael N. Alexander, Guillaume Dauphinais, J. Eli Bourassa,
#     Krishna K. Sabapathy, Nicolas C. Menicucci, and Ish Dhand (2021) *Fault-Tolerant Quantum Computation with Static Linear Optics*,
#     `PRX Quantum, Vol. 2, No. 4
#     <http://dx.doi.org/10.1103/PRXQuantum.2.040353>`__.
#
# .. [#ShorQEC1995]
#
#     Peter W. Shor (1995) *Scheme for reducing decoherence in quantum computer memory*,
#     `Physical Review A, Vol. 52, Iss. 4
#     <https://journals.aps.org/pra/abstract/10.1103/PhysRevA.52.R2493>`__.
#
# .. [#LatticeSurgeryRaussendorf2018]
#
#     Daniel Herr, Alexandru Paler, Simon J. Devitt and Franco Nori (2018) *Lattice Surgery on the Raussendorf Lattice*,
#     `IOP Publishing 3, 3
#     <https://arxiv.org/abs/1711.04921>`__.
#
# .. [#GoogleQEC2022]
#
#     Google Quantum AI (2022) *Quantum Error Correction in Quantum Computers*, `arXiv <https://arxiv.org/pdf/2207.06431.pdf>`__.
#
# .. [#FowlerSurfaceCode]
#
#     Austin G. Fowler, Matteo Mariantoni, John M. Martinis, Andrew N. Cleland (2012)
#     *Surface codes: Towards practical large-scale quantum computation*, `arXiv <https://arxiv.org/abs/1208.0928>`__.
