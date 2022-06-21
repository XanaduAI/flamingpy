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
""""Unit tests for GKP-specific functions in the GKP module."""

from datetime import datetime
import logging


import numpy as np
from numpy import sqrt, pi
from numpy.random import default_rng as rng
import pytest

from flamingpy.cv.gkp import integer_fractional, GKP_binner, Z_err, Z_err_cond


N = 50
now = datetime.now()
int_time = int(str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute))
logging.info("the following seed was used for random number generation: %i", int_time)

# Construct random numbers from an integer and fractional part.
alpha_vals = np.append(rng(int_time).random(5) * 5, np.sqrt(np.pi))


class TestGKPBinning:
    """Tests for GKP binning functions."""

    @pytest.mark.parametrize("alpha", alpha_vals)
    def test_integer_fractional(self, alpha):
        """Test that the integer and fractional part as obtained by
        integer_fractional matches that of constructed numbers."""
        integers = rng(int_time).integers(-N // 2, N // 2, N)
        fractions = (rng(int_time).random(N) - 0.5) * alpha
        numbers = integers * alpha + fractions
        int_part, frac_part = integer_fractional(numbers, alpha)
        assert np.all(int_part == integers)
        assert np.allclose(frac_part, fractions)

    def test_gkp_binner(self):
        """Tests that GKP_binner gives the integer part mod 2, and returns the
        fractional part if asked."""
        alpha = np.sqrt(np.pi)
        integers = rng(int_time).integers(-N // 2, N // 2, N)
        fractions = rng(int_time).random(N) * (alpha / 2)
        numbers = integers * alpha + fractions

        bits = integers % 2
        # int_part, frac_part = integer_fractional(numbers, alpha)
        assert np.all(GKP_binner(numbers) == bits)
        assert np.allclose(GKP_binner(numbers, return_fraction=True)[1], fractions)


# Even and odd homodyne outcomes mod sqrt(pi), and outcomes in the middle.
even_homs = np.array([2 * i * sqrt(pi) for i in range(-N // 2, N // 2)])
odd_homs = np.array([(2 * i + 1) * sqrt(pi) for i in range(-N // 2, N // 2)])
middle_homs = np.array([(2 * i + 1) * sqrt(pi) / 2 for i in range(-N // 2, N // 2)])
# Limit for summations.
lim = int(2 * N * np.sqrt(np.pi))

# Random low and high delta values
low_delta = rng(int_time).uniform(0.0001, 0.001, N)
high_delta = rng(int_time).uniform(10, 15, N)


class TestPhaseProbs:
    """Test the phase error proability functions."""

    def test_Z_err(self):
        """Ensure phase errors are 0 for low deltas and 0.5 for high deltas."""
        low_probs = Z_err(low_delta)
        high_probs = Z_err(high_delta)
        assert np.allclose(low_probs, 0)
        assert np.allclose(high_probs, 0.5)

    @pytest.mark.parametrize("use_hom_val", sorted([False, True]))
    def test_Z_err_cond(self, use_hom_val):
        """Test high-squeezing (low delta) regime."""
        for delta in (high_delta, low_delta):
            even_probs = Z_err_cond(delta, even_homs, var_num=lim, use_hom_val=use_hom_val)
            odd_probs = Z_err_cond(delta, odd_homs, var_num=lim, use_hom_val=use_hom_val)
            mid_probs = Z_err_cond(delta, middle_homs, var_num=lim, use_hom_val=use_hom_val)
            # Ensure that conditional phase error probabilities are close
            # to 0 for low delta values and 0.5 for high delta values.
            if np.array_equal(delta, high_delta):
                ref_prob = 0.5
            else:
                ref_prob = 0
            assert np.allclose(even_probs, ref_prob, rtol=0.001)
            assert np.allclose(odd_probs, ref_prob, rtol=0.001)
            assert np.allclose(mid_probs, ref_prob, rtol=0.001)

        # TODO: Test use_hom_val argument, changing summation limit,
        # 0-denominator behaviour
