FlamingPy Documentation
=======================

.. rst-class:: lead grey-text ml-2

:Release: |release|

.. raw:: html

    <style>
        .breadcrumb {
            display: none;
        }
        h1 {
            text-align: center;
            margin-bottom: 15px;
        }
        p.lead.grey-text {
            margin-bottom: 30px;
        }
        .footer-relations {
            border-top: 0px;
        }
    </style>

    <div class="container mt-2 mb-2">
       <div class="row container-fluid">
            <div class="col-lg-4 col-12 align-middle mb-2 text-center">
                <img src="_static/flamingpy_logo_light.svg" class="img-fluid" alt="Responsive image" style="width:100%; max-width: 300px;"></img>
            </div>
            <div class="col-lg-8 col-12 align-middle mb-2">
                <p class="lead grey-text">
                    FlamingPy (FP) is a cross-platform Python library with a variety of backends
                    for efficient simulations of error correction in fault-tolerant quantum computers.
                </p>
            </div>
        </div>
        <div class="row mt-3">

.. index-card::
    :name: Key Concepts
    :link: usage/getting_started.html
    :description: Learn how to simulate a fault-tolerant quantum computer

.. index-card::
    :name: Developing
    :link: development/contribution.html
    :description: How you can contribute to the FlamingPy project

.. index-card::
    :name: API
    :link: source/fp.html
    :description: Explore the FlamingPy API

.. raw:: html

        </div>
    </div>

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


Community
---------

You can start discussions and connect with our community members on our `Discussions Forum <https://discuss.pennylane.ai/c/flamingpy>`_.

Attribution for authors
-----------------------

FlamingPy is the work of `many contributors <https://github.com/XanaduAI/flamingpy/graphs/contributors>`__. If you are doing research using FlamingPy, please cite our paper below:

    Ilan Tzitrin, Takaya Matsuura, Rafael N. Alexander, Guillaume Dauphinais, J. Eli Bourassa, Krishna K. Sabapathy, Nicolas C. Menicucci, and Ish Dhand,
    Fault-Tolerant Quantum Computation with Static Linear Optics, PRX Quantum, Vol. 2, No. 4, 2021,
    |PRX_Quantum|

In addition to the authors above, the developers would like to thank Sanchit Bapat, Ashlesha Patil, Michael Vasmer, and Trevor Vincent for their contributions to the pre-release project.

License
-------

FlamingPy is **free** and **open source**, and released under the |Apache2|.


.. toctree::
   :maxdepth: 2
   :caption: Home
   :hidden:

   self
   install

.. toctree::
   :maxdepth: 2
   :caption: Background
   :hidden:

   quantum_error_correction

.. toctree::
   :maxdepth: 2
   :caption: Using FlamingPy
   :hidden:

   usage/getting_started
   usage/tutorials

.. toctree::
   :maxdepth: 2
   :caption: Development
   :hidden:

   development/guide_for_devs
   build_docs
   development/contribution

.. toctree::
   :maxdepth: 2
   :caption: Getting Help
   :hidden:

   faq
   help/support

.. toctree::
   :maxdepth: 2
   :caption: Python API
   :hidden:

   source/fp
   source/fp.codes
   source/fp.cv
   source/fp.decoders
   source/fp.noise
   source/fp.simulations
   source/fp.utils


.. |PRX_Quantum| raw:: html

   <a href="http://dx.doi.org/10.1103/PRXQuantum.2.040353" target="_blank">DOI:10.1103/prxquantum.2.040353</a>


.. |Apache2| raw:: html

   <a href="http://www.apache.org/licenses/LICENSE-2.0" target="_blank">Apache License, Version 2.0</a>