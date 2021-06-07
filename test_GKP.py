# Copyright 2020 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
""""Unit tests for GKP-specific functions in GKP.py."""
import pytest
import math
import numpy as np
from numpy import sqrt, pi
from numpy.random import default_rng as rng
from GKP import to_pi_string, integer_fractional, GKP_binner, Z_err, Z_err_cond

N = 50


def test_to_pi_string():
    """Test for the convenience function to_pi_string."""
    # Test +- sqrt(pi) and sqrt(pi)/2.
    assert to_pi_string(np.sqrt(np.pi)) == "$\\sqrt{\\pi}$"
    assert to_pi_string(-np.sqrt(np.pi)) == "$-\\sqrt{\\pi}$"
    assert to_pi_string(np.sqrt(np.pi) / 2) == "$\\sqrt{\\pi}/2$"
    assert to_pi_string(-np.sqrt(np.pi) / 2) == "$-\\sqrt{\\pi}/2$"
    # Test random odd integer multiples of sqrt(pi)/2.
    odd_int = 2 * rng().integers(-25, 25) - 1
    odd_int_label = "" if odd_int in (-1, 1) else odd_int
    assert to_pi_string(odd_int * np.sqrt(np.pi) / 2) == "${}\\sqrt{{\\pi}}/2$".format(
        odd_int_label
    )
    #  Test random even multiples of sqrt(pi).
    even_int = odd_int + 1
    assert to_pi_string(even_int * np.sqrt(np.pi)) == "${}\\sqrt{{\\pi}}$".format(
        even_int
    )
    # Check everything else converted into a str.
    rand_numb = rng().random()
    if not np.isclose(math.remainder(rand_numb, np.sqrt(np.pi) / 2), 0):
        assert to_pi_string(rand_numb) == str(rand_numb)


alpha = np.sqrt(np.pi)
integers = rng().integers(-25, 25, N)
fractions = rng().random(N) * (alpha / 2)
numbers = integers * alpha + fractions


class TestGKPBinning:
    """Tests GKP binning functions."""

    def test_integer_fractional(self):
        int_part, frac_part = integer_fractional(numbers, alpha)
        assert np.all(int_part == integers)
        assert np.allclose(frac_part, fractions)

    def test_gkp_binner(self):
        bits = integers % 2
        assert np.all(GKP_binner(numbers) == bits)
        assert np.allclose(GKP_binner(numbers, return_fraction=True)[1], fractions)


even_homs = np.array([2 * i * sqrt(pi) for i in range(-N // 2, N // 2)])
odd_homs = np.array([(2 * i + 1) * sqrt(pi) for i in range(-N // 2, N // 2)])
middle_homs = np.array([(2 * i + 1) * sqrt(pi) / 2 for i in range(-N // 2, N // 2)])
lim = int(2 * N * np.sqrt(np.pi))

low_delta = 0.0001 * (rng().random(N) + 1)
high_delta = np.full(N, 10)


class TestPhaseProbs:
    """Test the phase error proability functions."""

    def test_Z_err(self):
        low_probs = Z_err(low_delta)
        high_probs = Z_err(high_delta)
        assert np.allclose(low_probs, 0)
        assert np.allclose(high_probs, 0.5)

    @pytest.mark.parametrize("use_hom_val", [False, True])
    def test_Z_err_cond(self, use_hom_val):
        # Test high-squeezing (low delta) regime.
        for delta in (high_delta, low_delta):
            even_probs = Z_err_cond(
                delta, even_homs, var_num=lim, use_hom_val=use_hom_val
            )
            odd_probs = Z_err_cond(
                delta, odd_homs, var_num=lim, use_hom_val=use_hom_val
            )
            mid_probs = Z_err_cond(
                delta, middle_homs, var_num=lim, use_hom_val=use_hom_val
            )
            print(delta)
            if np.array_equal(delta, high_delta):
                ref_prob = 0.5
            else:
                ref_prob = 0
            print(ref_prob, even_probs)
            assert np.allclose(even_probs, ref_prob, rtol=0.001)
            assert np.allclose(odd_probs, ref_prob, rtol=0.001)
            assert np.allclose(mid_probs, ref_prob, rtol=0.001)

        # TODO: Test use_hom_val argument, changing summation limit,
        # 0-denominator behaviour
