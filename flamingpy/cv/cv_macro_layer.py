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
"""Functions for reducing a macronode lattice to a canonical lattice."""

# pylint: disable=protected-access

import numpy as np
from scipy.linalg import block_diag

from flamingpy.cv.ops import CVLayer, SCZ_apply
from flamingpy.cv.gkp import GKP_binner, Z_err_cond
from thewalrus.symplectic import expand, beam_splitter


