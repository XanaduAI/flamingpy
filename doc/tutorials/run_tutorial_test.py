r"""
.. _run_tutorial_test:
Basic tutorial: qubit rotation
==============================

*Author: PennyLane dev team. Last updated: 19 Jan 2021.*
To see how PennyLane allows the easy construction and optimization of quantum functions, let's
consider the simple case of **qubit rotation** the PennyLane version of the 'Hello, world!'
example.
The task at hand is to optimize two rotation gates in order to flip a single
qubit from state :math:`\left|0\right\rangle` to state :math:`\left|1\right\rangle`.
"""

##############################################################################
#
# Let's see how we can easily implement and optimize this circuit using PennyLane.
#
# Importing PennyLane and NumPy
# -----------------------------
#
# The first thing we need to do is import PennyLane, as well as the wrapped version
# of NumPy provided by PennyLane.

import flamingpy as fp
import numpy as np
import matplotlib.pyplot as plt


##############################################################################
# .. important::
#
#     When constructing a hybrid quantum/classical computational model with PennyLane,
#     it is important to **always import NumPy from PennyLane**, not the standard NumPy!
#
#     By importing the wrapped version of NumPy provided by PennyLane, you can combine
#     the power of NumPy with PennyLane:
#
#     * continue to use the classical NumPy functions and arrays you know and love
#     * combine quantum functions (evaluated on quantum hardware/simulators) and
#       classical functions (provided by NumPy)
#     * allow PennyLane to automatically calculate gradients of both classical and
#       quantum functions

x, y = np.random.rand(10), np.random.rand(10)
plt.scatter(x,y)