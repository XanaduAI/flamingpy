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
"""Benchmark Lemon's max_weight_matching against NetworkX."""
import time

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from flamingpy.decoders.mwpm.lemon import max_weight_matching


time_lemon = []
time_nx = []
matrix_size = []

# Change this for more data points
for num_nodes in range(50, 150, 50):

    matrix_size.append(num_nodes)
    b = np.random.randint(0, 2001, size=(num_nodes, num_nodes))
    b_symm = (b + b.T) * 0.5
    for i in range(num_nodes):
        b_symm[i, i] = 0

    G_match = nx.from_numpy_array(b_symm)

    start = time.time()
    lemon_matching = max_weight_matching(G_match, weight="weight")
    end = time.time()
    time_lemon.append(end - start)

    start = time.time()
    nx_matching = nx.max_weight_matching(G_match, weight="weight")
    end = time.time()
    time_nx.append(end - start)

    print(
        "Matrix size: {}; t_lemon: {}; t_nx: {}".format(
            matrix_size[-1], time_lemon[-1], time_nx[-1]
        )
    )

plt.plot(matrix_size, time_lemon, label="lemon")
plt.plot(matrix_size, time_nx, label="NetworkX")
plt.xlabel("Nodes")
plt.ylabel("Time")
plt.yscale("log")
plt.legend(loc="best")
if __name__ == "__main__":
    plt.savefig("lemon_benchmark.pdf")
