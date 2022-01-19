from ft_stack import viz
from ft_stack import decoder as dec
from ft_stack.RHG import RHGCode, alternating_polarity
from ft_stack.graphstates import CVGraph
import matplotlib.pyplot as plt


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
    "delta": delta, }

# Drawing options
dw = {
    "show_nodes": True,
    "color_nodes": "state",
    "label": "bit_val",
    "legend": True,
    "title": True,
    "display_axes": True,
    "label_edges": True,
    "label_cubes": True,
    "label_boundary": False,
}

# Apply noise
CVRHG.apply_noise(cv_noise)
# Measure syndrome
CVRHG.measure_hom("p", RHG_code.syndrome_inds)

# Manual decoding to plot intermediate results.
dec.CV_decoder(RHG_code, translator=dec.GKP_binner)
G_dec, G_match = dec.build_dec_and_match_graphs(RHG_code, weight_options)
matching = G_match.min_weight_perfect_matching()
viz.draw_dec_graph(G_dec, dw.get("label_edges"))
ax = viz.syndrome_plot(
    RHG_code, G_dec, index_dict=RHG_code._decoder_mapping, drawing_opts=dw
)
viz.draw_matching_on_syndrome_plot(ax, matching, G_dec, G_match, dw.get("label_edges"))
# This function requires a network graph object. Most backends implement
# the to_nx() method to perform the conversion if needed.
viz.draw_dec_graph(G_match.graph, title="Matching graph")
plt.show()

# Automatic decoding
c = dec.correct(
    code=RHG_code,
    decoder=decoder,
    weight_options=weight_options,
    sanity_check=True,
)
print(f"Success: {c}")
