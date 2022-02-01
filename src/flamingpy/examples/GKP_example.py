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

import numpy as np
from flamingpy import GKP, viz


alpha = np.sqrt(np.pi)
xs = np.arange(-10, 10, 0.01)

ns, fs = GKP.integer_fractional(xs, alpha)
viz.plot_integer_fractional(xs, ns, fs, alpha)

bit_values = GKP.GKP_binner(xs)
viz.plot_GKP_bins(xs, bit_values, alpha)

delta = 0.1
error_hom_val = GKP.Z_err_cond([delta] * len(xs), xs, use_hom_val=True)
error_no_hom_val = GKP.Z_err_cond([delta] * len(xs), xs)

viz.plot_Z_err_cond(xs, error_hom_val, alpha, True)
viz.plot_Z_err_cond(xs, error_no_hom_val, alpha, False)
