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
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import matplotlib.pyplot as plt


def to_pi_string(x, tex=1):
    """Convert x, a multiple of sqrt(pi)/2, to a string."""
    remainder = math.remainder(x, np.sqrt(np.pi) / 2)
    if not round(remainder):
        integer = round(x / (np.sqrt(np.pi)/2))
        pref = int(integer * ((1 - integer % 2) / 2 + integer % 2))
        x_str = (not bool(round(x))) * '0' + bool(round(x)) * (
            bool(tex) * '$' + (not bool(1 + pref)) * '-' +
            bool(1 - abs(pref)) * str(pref) + r'\sqrt{\pi}' + (integer % 2) * '/2' +
            bool(tex) * '$'
            )
        return x_str
    return str(x)

def basic_translate(outcomes):
    """Naively translate CV outcomes to bit values.

    The function treats values in (-sqrt(pi)/2, sqrt(pi)/2) as 0
    and those in (sqrt(pi)/2, 3sqrt(pi)/2) as 1. Values on the
    boundary are assigned a random bit value. The rest of the bins
    are defined periodically.
    Args:
        outcomes (array): the values of a p-homodyne measurement
    Retruns:
        array: the corresponding bit values.
    """
    # Bin width
    alpha = np.sqrt(np.pi)
    n = len(outcomes)
    bit_values = np.zeros(n, dtype=int)
    for i in range(n):
        div = np.divmod(outcomes[i], alpha)
        if div[1] < alpha/2:
            bit_values[i] = div[0] % 2
        elif div[1] > alpha / 2:
            bit_values[i] = (div[0]+1) % 2
        else:
            bit_values[i] = np.random.randint(2)
    return bit_values


def Z_err(var, var_num=5):
    """Return the probability of Z errors for a list of variances.

    Args:
        var (array): array of lattice p variances
        var_num (float): number of variances away from the origin we
            include in the integral
        translator (function): the CV-to-bit translator; the basic
            binning function by default
    Returns:
        array: probability of Z (phase flip) errors for each variance.
    """
    # Find largest bin number by finding largest variance, multiplying by
    # var_num, then rounding up to nearest integer of form 4n+1, which are
    # the left boundaries of the 0 bins mod sqrt(pi)
    n_max = int(np.ceil(var_num*np.amax(var))//2*4 + 1)
    # error = 1 - integral over the 0 bins
    # Initiate a list with length same as var
    error = np.ones(len(var))
    # Integral over 0 bins that fell within var_num*var_max away from
    # origin
    for i in range(-n_max, n_max, 4):
        error -= 0.5*(erf((i+2)*np.sqrt(np.pi)/(2*var))
                      - erf(i*np.sqrt(np.pi)/(2*var)))
    return error

def Z_err_cond(var, hom_val, var_num=5):
    """Return the conditional phase error probability for lattice nodes.

    Return the phase error probability for a list of variances var
    given homodyne outcomes hom_val, with var_num used to determine 
    the number of terms to keep in the summation in the formula.

    Args:
        var (array): the lattice p variances
        hom_val (array): the p-homodyne outcomes
        var_num (float): number of variances away from the origin we
            include in the integral
    Returns:
        array: probability of Z (phase flip) errors for each variance,
            contioned on the homodyne outcomes.
    """
    # Find largest bin number by finding largest variance, multiplying by
    # var_num, then rounding up to nearest integer of form 4n+1, which are
    # the left boundaries of the 0 bins mod sqrt(pi)
    n_max = int(np.ceil(var_num * np.amax(var)) // 2 * 4 + 1)
    # Initiate a list with length same as var
    # TODO replace ex with normal pdf?
    ex = lambda z, n: np.exp(-(z - n * np.sqrt(np.pi)) ** 2 / var)
    error = np.zeros(len(var))
    alpha = np.sqrt(np.pi)
    # Take modulus to obtain a new range of 0 to sqrt(pi).
    mod_val = np.mod(hom_val, alpha)
    # If mod_val is less than sqrt(pi)/2, keep it as is. If greater, convert
    # it to sqrt(pi) - mod_val.
    z = mod_val + (((np.sign(alpha/2 - mod_val) - 1)) / 2) * alpha
    numerator = np.sum([ex(z, 2*i+1) for i in range(-n_max, n_max)], 0)
    denominator = np.sum([ex(z, i) for i in range(-n_max, n_max)], 0)
    error = numerator / denominator
    return error.round(5)


if __name__ == '__main__':
    pass