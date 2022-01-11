import time

import matplotlib.pyplot as plt

from ft_stack import decoder as dec
from ft_stack import matching as mt
from ft_stack.graphstates import CVGraph
from ft_stack.RHG import RHGCode, alternating_polarity

# How many simulations to do for each algorithm.
num_trials = 10

# DV (outer) code
distance = 3
boundaries = "periodic"
RHG_code = RHGCode(
    distance=distance, boundaries=boundaries, polarity=alternating_polarity
)
RHG_lattice = RHG_code.graph
# CV (inner) code/state
p_swap = 0
CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)

# Noise model
delta = 0.1
cv_noise = {"noise": "grn", "delta": delta, "sampling_order": "initial"}

# Decoding options
decoder = {"inner": "basic", "outer": "MWPM"}
weight_options = {
    "method": "blueprint",
    "integer": True,
    "multiplier": 100,
    "delta": delta,
}

times = {
    "networkx": [],
    "lemon": [],
    "retworkx": [],
}

matching_graph = {
    "networkx": mt.NxMatchingGraph,
    "lemon": mt.LemonMatchingGraph,
    "retworkx": mt.RxMatchingGraph,
}


for i in range(num_trials):
    for alg in ["networkx", "lemon", "retworkx"]:
        # Apply noise
        CVRHG.apply_noise(cv_noise)
        # Measure syndrome
        CVRHG.measure_hom("p", RHG_code.syndrome_inds)

        # Manual decoding to plot intermediate results.
        dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
        G_dec, G_match = dec.build_dec_and_match_graphs(
            RHG_code, weight_options, matching_backend=matching_graph[alg]
        )
        before = time.time()
        matching = G_match.min_weight_perfect_matching()
        after = time.time()
        times[alg].append(after - before)

plt.figure()
# bins = np.logspace(-3, 0, 30)
for alg in ["networkx", "lemon", "retworkx"]:
    plt.hist(times[alg], bins=10, label=alg)
plt.legend()
plt.xscale("log")
plt.xlabel("Times [seconds]")
plt.ylabel("Count")
plt.title(f"Matching for code distance {distance}")
plt.savefig(f"benchmark_matching_distance_{distance}.pdf")
