from ft_stack import viz
from ft_stack.graphstates import EGraph, CVGraph
import matplotlib.pyplot as plt

# Bell state EGraph
edge = [(0, 0, 0), (1, 1, 1)]
dims = (1, 1, 1)
bell_state = EGraph(dims=dims)
bell_state.add_edge(*edge, color="MidnightBlue")
# Plot the bell state
viz.draw_EGraph(
    bell_state, color_edges="MidnightBlue", color_nodes="magenta", label="index"
)
plt.show()
bell_state.adj_generator(sparse=True)
print("Adjacency matrix: \n", bell_state.adj_mat, "\n")

CVbell = CVGraph(bell_state, p_swap=0.5)
# Noise model for CVGraph
model = {"noise": "grn", "delta": 1, "sampling_order": "final"}
CVbell.apply_noise(model)
CVbell.measure_hom("p", [0])
CVbell.measure_hom("q", [1])
CVbell.eval_Z_probs(cond=False)

# Some parameters to make the graph of CVbell.
CV_graph_params = {
    "color_nodes": "state",
    "legend": True,
    "title": True,
    "state_colors": {state: None for state in CVbell._states},
}

viz.draw_EGraph(CVbell.egraph, label="hom_val_p", **CV_graph_params)
viz.draw_EGraph(CVbell.egraph, label="hom_val_q", **CV_graph_params)
viz.draw_EGraph(CVbell.egraph, label="p_phase", **CV_graph_params)
plt.show()

print("\nNodes :", bell_state.nodes.data())
print("Edges :", bell_state.edges.data())
print("p indices: ", CVbell.p_inds, "\n")
print("GKP indices: ", CVbell.GKP_inds, "\n")
print("\nSymplectic CZ matrix: \n", CVbell.SCZ(), "\n")

viz.plot_SCZ_mat_heat_map(CVbell.SCZ())
