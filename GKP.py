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
    """Convert x, a multiple of sqrt(pi)/2, to a pretty string.

    If x is not a multiple of sqrt(pi)/2, return the unmodified string
    of x. If tex is True, add LaTeX $ signs.
    """
    remainder = math.remainder(x, np.sqrt(np.pi) / 2)
    if np.isclose(remainder, 0):
        integer = round(x / (np.sqrt(np.pi) / 2))
        pref = int(integer * ((1 - integer % 2) / 2 + integer % 2))
        x_str = (not bool(round(x))) * "0" + bool(round(x)) * (
            bool(tex) * "$"
            + (not bool(1 + pref)) * "-"
            + bool(1 - abs(pref)) * str(pref)
            + r"\sqrt{\pi}"
            + (integer % 2) * "/2"
            + bool(tex) * "$"
        )
        return x_str
    return str(x)


vec_remainder = np.vectorize(math.remainder)


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
    f = vec_remainder(x, alpha)
    # The integer part. astype(int) used here to prevent unwanted
    # behaviour during binning (e.g. -6.01 % 2 = 2, while -6 % 2 =0)
    n = ((x - f) / alpha).astype(int)
    if draw:
        xmin, xmax = alpha * (x[0] // alpha), alpha * (x[-1] // alpha) + alpha
        newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
        newxlabels = [to_pi_string(tick) for tick in newxticks]
        plt.plot(x, n, ",")
        plt.title("Integer Part")
        plt.xticks(newxticks, newxlabels, fontsize="small")
        plt.show()

        plt.title("Fractional Part")
        plt.plot(x, f, ",")
        newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
        newylabels = ["{:.3f}".format(tick) for tick in newyticks[1:-1]]
        newylabels = [to_pi_string(-alpha / 2)] + newylabels + [to_pi_string(alpha / 2)]
        plt.xticks(newxticks, newxlabels, fontsize="small")
        plt.yticks(newyticks, newylabels)
        plt.show()
    return n, f


def GKP_binner(outcomes, return_fraction=False, draw=False):
    """Naively translate CV outcomes to bit values.

    The function treats values in (-sqrt(pi)/2, sqrt(pi)/2) as 0
    and those in (sqrt(pi)/2, 3sqrt(pi)/2) as 1. Values on the
    boundary are binned to 0. The rest of the bins
    are defined periodically.

    Args:
        outcomes (array): the values of a p-homodyne measurement.
        return_fraction (bool): return the fractional part of the
            outcome as well, if desired.
        draw (bool): if True, sketch binned values over outcomes.

    Retruns:
        array: the corresponding bit values.
    """
    # Bin width
    alpha = np.sqrt(np.pi)
    # np.divmod implementation also possible
    int_frac = integer_fractional(outcomes, alpha)
    # CAUTION: % has weird behaviour that seems to only show up for
    # large (n ~ 100) multiples of sqrt(pi). Check!
    bit_values = int_frac[0] % 2
    if draw:
        xmin, xmax = alpha * (x[0] // alpha), alpha * (x[-1] // alpha) + alpha
        newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
        newxlabels = [to_pi_string(tick) for tick in newxticks]
        plt.plot(outcomes, bit_values, ",")
        plt.title("Binned values")
        plt.xticks(newxticks, newxlabels, fontsize="small")
        plt.show()
    if return_fraction:
        return bit_values, int_frac[1]
    return bit_values


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
        error -= 0.5 * (
            erf((i + 2) * np.sqrt(np.pi) / (2 * var))
            - erf(i * np.sqrt(np.pi) / (2 * var))
        )
    return error


def Z_err_cond(var, hom_val, var_num=10, replace_undefined=0, use_hom_val=False):
    """Return the conditional phase error probability for lattice nodes.

    Return the phase error probability for a list of variances var
    given homodyne outcomes hom_val, with var_num used to determine
    the number of terms to keep in the summation in the formula.
    By default, the fractional part of hom_val is used in the summation;
    if use_hom_val is True, use the entire homodyne value.

    Args:
        var (array): the lattice p variances.
        hom_val (array): the p-homodyne outcomes.
        var_num (float): number of variances away from the origin we
            include in the integral.
        replace_undefined (float): how to handle 0 denominators.
            If 'bin_location', return poor-man's probability that ranges
            from 0 in the centre of a bin to 0.5 halfway between bins.
            Otherwise, sets it to the replace_undefined value (0 by
            default).
        use_hom_val (bool): if True, use the entire homodyne value
            hom_val in the expression; otherwise use the fractional
            part.
    Returns:
        array: probability of Z (phase flip) errors for each variance,
            contioned on the homodyne outcomes.
    """
    # TODO: Make the following line smarter.
    n_max = var_num
    # Initiate a list with length same as var
    # TODO replace ex with normal pdf?
    ex = lambda z, n: np.exp(-((z - n * np.sqrt(np.pi)) ** 2) / var)
    error = np.zeros(np.shape(var))
    bit, frac = GKP_binner(hom_val, return_fraction=True)
    factor = 1 - bit
    val = hom_val if use_hom_val else frac
    numerator = np.sum([ex(val, 2 * i + factor) for i in range(-n_max, n_max)], 0)
    denominator = np.sum([ex(val, i) for i in range(-n_max, n_max)], 0)
    # Dealing with 0 denonimators
    where_0 = np.where(denominator == 0)[0]
    the_rest = np.delete(np.arange(np.size(var)), where_0)
    error = np.empty(np.size(var))

    # For 0 denominator, populate error according to replace_undefined
    if replace_undefined == 'bin_location':
        zero_dem_result = np.abs(val) / np.sqrt(np.pi)
    else:
        zero_dem_result = replace_undefined
    error[where_0] = np.full(np.len(where_0), zero_dem_result)

    if np.size(var) == 1:
        numerator = np.array([numerator])
        denominator = np.array([denominator])
    error[the_rest] = numerator[the_rest] / denominator[the_rest]
    if np.size(var) == 1:
        return error[0]
    return error


if __name__ == "__main__":
    alpha = np.sqrt(np.pi)
    x = np.arange(-10, 10, 0.01)
    integer_fractional(x, alpha, draw=True)
    GKP_binner(x, draw=True)
