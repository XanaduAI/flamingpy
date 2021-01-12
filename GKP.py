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
"""GKP-state-specific functions."""
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf


def to_pi_string(x, tex=True):
    """Convert x, a multiple of sqrt(pi)/2, to a string.

    If tex is True, convert to LaTeX.
    """
    remainder = math.remainder(x, np.sqrt(np.pi) / 2)
    if not round(remainder):
        integer = round(x / (np.sqrt(np.pi) / 2))
        pref = int(integer * ((1 - integer % 2) / 2 + integer % 2))
        x_str = (
            not bool(round(x))) * '0' + bool(round(x)) * (
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
        if div[1] < alpha / 2:
            bit_values[i] = div[0] % 2
        elif div[1] > alpha / 2:
            bit_values[i] = (div[0] + 1) % 2
        else:
            bit_values[i] = np.random.randint(2)
    return bit_values


def integer_fractional(x, alpha, draw=False):
    """Obtain the integer and fractional part of x with respect to alpha.

    Any real number x can be expressed as n * alpha + f, where n is an
    integer, alpha is any real number, and f is a real number such that
    |f| <= alpha / 2. This function returns n and f.

    Args:
        x (float): a real number
        alpha (float): alpha from above
        draw (bool): if True, plot the fractional and integer parts
            of x.
    """
    # The fractional part.
    f = np.vectorize(math.remainder)(x, alpha)
    # The integer part.
    n = (x - f) / alpha
    if draw:
        xmin, xmax = alpha * (x[0] // alpha), alpha * (x[-1] // alpha) + alpha
        newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
        newxlabels = [to_pi_string(tick) for tick in newxticks]
        plt.plot(x, n, ',')
        plt.title('Integer Part')
        plt.xticks(newxticks, newxlabels, fontsize='small')
        plt.show()

        plt.title('Fractional Part')
        plt.plot(x, f, ',')
        newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
        newylabels = ['{:.3f}'.format(tick) for tick in newyticks[1:-1]]
        newylabels = [to_pi_string(-alpha / 2)] + newylabels + [to_pi_string(alpha / 2)]
        plt.xticks(newxticks, newxlabels, fontsize='small')
        plt.yticks(newyticks, newylabels)
        plt.show()
    return n, f


def Z_err(var, var_num=5):
    """Return the probability of Z errors for a list of variances.

    Args:
        var (array): array of lattice p variances
        var_num (float): number of variances away from the origin we
            include in the integral
    Returns:
        array: probability of Z (phase flip) errors for each variance.
    """
    # Find largest bin number by finding largest variance, multiplying by
    # var_num, then rounding up to nearest integer of form 4n+1, which are
    # the left boundaries of the 0 bins mod sqrt(pi)
    n_max = int(np.ceil(var_num * np.amax(var)) // 2 * 4 + 1)
    # error = 1 - integral over the 0 bins
    # Initiate a list with length same as var
    error = np.ones(len(var))
    # Integral over 0 bins that fell within var_num*var_max away from
    # origin
    for i in range(-n_max, n_max, 4):
        error -= 0.5 * (erf((i + 2) * np.sqrt(np.pi) / (2 * var))
                      - erf(i * np.sqrt(np.pi) / (2 * var)))
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
    # TODO: Make the following line smarter.
    n_max = var_num
    # Initiate a list with length same as var
    # TODO replace ex with normal pdf?
    ex = lambda z, n: np.exp(-(z - n * np.sqrt(np.pi)) ** 2 / var)
    error = np.zeros(len(var))
    z = integer_fractional(hom_val, np.sqrt(np.pi))[1]
    numerator = np.sum([ex(z, 2 * i + 1) for i in range(-n_max, n_max)], 0)
    denominator = np.sum([ex(z, i) for i in range(-n_max, n_max)], 0)
    error = numerator / denominator
    return error.round(5)


if __name__ == '__main__':
    alpha = np.sqrt(np.pi)
    x = np.arange(-5, 5, 0.01)
    integer_fractional(x, alpha, draw=True)
