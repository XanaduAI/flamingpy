# Copyright 2022 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""This subpackage contains the modules that make up the Flamingpy decoders. This includes the
MWPM and Unionfind decoders.

Available modules
-----------------

.. currentmodule:: flamingpy.decoders
.. autosummary::
    :toctree: api

    ~decoder
    ~mwpm
    ~unionfind

Module members
--------------

.. currentmodule:: flamingpy.decoders
.. autosummary::
    :toctree: api

    correct
    uf_decoder
    mwpm_decoder

Class Inheritance Diagram
--------------------------
.. inheritance-diagram:: flamingpy.decoders.mwpm.LemonMatchingGraph
     flamingpy.decoders.mwpm.RxMatchingGraph
   :parts: 1
"""
from .decoder import correct
from .unionfind.algos import uf_decoder
from .mwpm.algos import mwpm_decoder
