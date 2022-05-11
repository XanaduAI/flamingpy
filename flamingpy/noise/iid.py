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
"""An implementation of IID Pauli noise."""

from numpy.random import default_rng
import numpy as np


class IidNoise:
    """Noise sampler for independent and identically distributed Z errors on
    the qubits of a cluster state.

    Args:
        code (SurfaceCode): the code on which to apply the noise.
        error_probability (float): the probability of a Z error.
    """

    def __init__(self, length, error_probability):
        if error_probability < 0.0 or error_probability > 1.0:
            raise ValueError("Probability is not between 0 and 1.")
        self.length = length
        self.error_probability = error_probability

    def sample(self, rng=default_rng()):
        """Generate a random error.

        Args:
            rng (numpy.random.Generator, optional): a random number generator
                following the NumPy API. It can be seeded for reproducibility.
                By default, numpy.random.default_rng is used without a fixed seed.

        Returns (numpy.ndarray):
            A bit string where each 1 represents
            a Z error on the corresponding qubit.
        
        """
        return np.array(rng.random(self.length) < self.error_probability, dtype=int)
