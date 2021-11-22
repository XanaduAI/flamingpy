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
""" 
Helper functions to draw various graphs 
and generate plots using Matplotlib.
"""
import matplotlib.pyplot as plt
import numpy as np
import ft_stack.GKP as GKP


def plot_integer_fractional(xs, ns, fs, alpha):
    xmin, xmax = alpha * (xs[0] // alpha), alpha * (xs[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(xs, ns, ",")
    plt.title("Integer Part", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.show()

    plt.title("Fractional Part", fontsize="medium")
    plt.plot(xs, fs, ",")
    newyticks = np.linspace(-alpha / 2, alpha / 2, num=7)
    newylabels = ["{:.3f}".format(tick) for tick in newyticks[1:-1]]
    newylabels = (
        [GKP.to_pi_string(-alpha / 2)] + newylabels + [GKP.to_pi_string(alpha / 2)]
    )
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks(newyticks, newylabels)
    plt.show()


def plot_GKP_bins(outcomes, bit_values, alpha):
    xmin, xmax = alpha * (outcomes[0] // alpha), alpha * (outcomes[-1] // alpha) + alpha
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(outcomes, bit_values, ",")
    plt.title("Binned values", fontsize="medium")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.yticks([0, 1], [0, 1])
    plt.show()


def plt_Z_err_cond(hom_val, error, alpha, use_hom_val):
    _, frac = GKP.GKP_binner(hom_val, return_fraction=True)
    val = hom_val if use_hom_val else frac
    xmin, xmax = alpha * (hom_val[0] // alpha), alpha * (hom_val[-1] // alpha) + alpha
    print(xmin, xmax, min(val), max(val))
    newxticks = np.linspace(xmin, xmax, int((xmax - xmin) // alpha) + 1)
    newxlabels = [GKP.to_pi_string(tick) for tick in newxticks]
    plt.plot(val, error, ",")
    addendum = "Full homodyne value" if use_hom_val else "Central peak"
    plt.title("Conditional phase probabilities: " + addendum, fontsize="small")
    plt.xticks(newxticks, newxlabels, fontsize="small")
    plt.show()
