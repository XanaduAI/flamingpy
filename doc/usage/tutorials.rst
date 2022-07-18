 .. role:: html(raw)
   :format: html

Tutorials
=========

..
    To add a tutorial, use the ``gallery-item`` directive from the ``xanadu-sphinx-theme``
    Note that the ``description`` option can be a link to a document. Also,
    thumbnails will be created out of figures and stored in ``tutorials/_out/images/thumb/``
    with the same name of the tutorial prepended with ``sphx_glr_``.
    Therefore, consider ``tutorials/_out`` as a "built" directory.

    **Example**

    .. gallery-item::
        :tooltip: This tutorial is directed at people who are new to FlamingPy.
        :figure: tutorials/_out/images/thumb/sphx_glr_run_intro_tutorial.png
        :description: :doc:`../tutorials/_out/run_intro_tutorial`

.. raw:: html

    <link href="https://cdnjs.cloudflare.com/ajax/libs/mdbootstrap/4.8.10/css/mdb.min.css" rel="stylesheet">

:html:`<div class="gallery-grid row">`

.. gallery-item::
    :tooltip: In this tutorial we will tell you a bit about graph states, and show you how to define and visualize them using FlamingPy.
    :figure: _static/graph_states_thumbnail.png
    :description: :doc:`../tutorials/_out/run_graph_states`

.. gallery-item::
    :tooltip: In this tutorial we will go over the minimal steps required to run through one round of quantum error correction: encoding, decoding, and recovery.
    :figure: _static/ec_thumbnail.png
    :description: :doc:`../tutorials/_out/run_error_correction`

.. gallery-item::
    :tooltip: Fault-tolerant measurement-based quantum computation with continuous-variable cluster states
    :figure: tutorials/ft_mbqc/mbqc_blueprint.png
    :description: :doc:`../tutorials/_out/run_ft_mbqc`

:html:`</div></div><div style='clear:both'>`
