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
"""The qubit codes module.

Abstract code modules
---------------------

.. currentmodule:: flamingpy.codes
.. autosummary::
    :recursive:
    :toctree: api

    graphs
    stabilizer

Specific code implementations
-----------------------------

.. currentmodule:: flamingpy.codes
.. autosummary::
    :recursive:
    :toctree: api

    surface_code
"""
from .stabilizer import Stabilizer
from .surface_code import alternating_polarity, RHG_graph, SurfaceCode
