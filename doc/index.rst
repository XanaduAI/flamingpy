FlamingPy Documentation
=======================

.. rst-class:: lead grey-text ml-2

:Release: |release|

.. raw:: html

    <style>
        #right-column.card {
            box-shadow: none!important;
        }
        #right-column.card:hover {
            box-shadow: none!important;
        }
        .breadcrumb {
            display: none;
        }
        h1 {
            text-align: center;
            margin-bottom: 15px;
        }
        .footer-relations {
            border-top: 0px;
        }
    </style>
    <div class="row container-fluid">
      <div class="col-lg-4 col-12 align-middle mb-2 text-center">
          <img src="_static/laptop.png" class="img-fluid" alt="Responsive image" style="width:100%; max-width: 300px;"></img>
      </div>
      <div class="col-lg-8 col-12 align-middle mb-2">
        <p class="lead grey-text">
            FlamingPy is a cross-platform Python library with a variety of backends 
            for efficient simulations of error correction in fault-tolerant quantum computers.
        </p>
      </div>
    </div>
    <div style='clear:both'></div>
    <div class="container mt-2 mb-2">
        <div class="row mt-3">
            <div class="col-lg-4 mb-2 adlign-items-stretch">
                <a href="usage/getting_started.html">
                    <div class="card rounded-lg" style="height:100%;">
                        <div class="d-flex">
                            <div>
                                <h3 class="card-title pl-3 mt-4">
                                Using SF
                                </h3>
                                <p class="mb-3 grey-text px-3">
                                    Learn how to simulate a fault-tolerant quantum computer <i class="fas fa-angle-double-right"></i>
                                </p>
                            </div>
                        </div>
                    </div>
                </a>
            </div>
            <div class="col-lg-4 mb-2 align-items-stretch">
                <a href="development/contribution.html">
                <div class="card rounded-lg" style="height:100%;">
                    <div class="d-flex">
                        <div>
                            <h3 class="card-title pl-3 mt-4">
                            Developing
                            </h3>
                            <p class="mb-3 grey-text px-3">How you can contribute to the FlamingPy project <i class="fas fa-angle-double-right"></i></p>
                        </div>
                    </div>
                </div>
            </a>
            </div>
            <div class="col-lg-4 mb-2 align-items-stretch">
                <a href="source/fp.html">
                <div class="card rounded-lg" style="height:100%;">
                    <div class="d-flex">
                        <div>
                            <h3 class="card-title pl-3 mt-4">
                            API
                            </h3>
                            <p class="mb-3 grey-text px-3">Explore the FlamingPy API <i class="fas fa-angle-double-right"></i></p>
                        </div>
                    </div>
                </div>
            </a>
            </div>
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


Attribution for authors
-----------------------

FlamingPy is the work of `many contributors <https://github.com/XanaduAI/ft-stack/graphs/contributors>`__. The developers would like to specifically thank Michael Vasmer and Ashlesha Patil for their contributions to the pre-release project. 

If you are doing research using FlamingPy, please cite our paper below:

    Ilan Tzitrin, Takaya Matsuura, Rafael N. Alexander, Guillaume Dauphinais, J. Eli Bourassa, Krishna K. Sabapathy, Nicolas C. Menicucci, and Ish Dhand,
    Fault-Tolerant Quantum Computation with Static Linear Optics, PRX Quantum, Vol. 2, No. 4, 2021, 
    `DOI:10.1103/prxquantum.2.040353 <http://dx.doi.org/10.1103/PRXQuantum.2.040353>`__ 


License
-------

FlamingPy is **free** and **open source**, and released under the `Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`__.


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
