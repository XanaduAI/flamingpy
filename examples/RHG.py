from ft_stack.RHG import RHGCode, alternating_polarity
from ft_stack import viz
import matplotlib.pyplot as plt

# Simple tests
# Instantiate an RHG latice of a certian distance, with certain
# boundaries. Draw the EGraph.
d = 2
# boundaries = "finite"
# boundaries = "primal"
boundaries = "periodic"
# boundaries = ["primal", "dual", "periodic"]
# For iterating through boundaries
# for boundaries in it.product(['primal', 'dual', 'periodic'], repeat=3):

RHG = RHGCode(d, boundaries=boundaries, polarity=alternating_polarity)
RHG_lattice = RHG.graph
# Check maronode lattice
# RHG_lattice = RHG_graph(d, boundaries=boundaries, macronodes=True)
ax = viz.draw_code_lattice(RHG_lattice)
plt.show()

# # Check edges between boundaries for periodic boundary conditions.
# all_boundaries = []
# for plane in ("x", "y", "z"):
#     for i in (0, 2 * d - 1):
#         all_boundaries += RHG.graph.slice_coords(plane, i)
#     # For macronode lattice, change to (-0.1, 2 * d - 0.9)
# RHG_subgraph = RHG_lattice.subgraph(all_boundaries)
# RHG_subgraph.draw()

# # Check stabilizer coordinates
# syndrome = RHG.stabilizers
# print("Number of six-body stabilizers :", len(syndrome))
# for i in range(len(syndrome)):
#     cube = syndrome[i]
#     color = np.random.rand(3)
#     for point in cube.egraph:
#         x, z, y = point
#         ax.scatter(x, z, y, color=color, s=200)
# ax.set_title(str(boundaries).capitalize() + " boundaries")

# # Check sampling
# delta = 0.001
# # Percent p-squeezed states.
# p_swap = 0
# CVRHG = CVGraph(RHG_lattice, p_swap=p_swap)
# for sampling_order in ["initial", "final", "two-step"]:
#     model = {"noise": "grn", "delta": delta, "sampling_order": sampling_order}
#     CVRHG.apply_noise(model)
#     CVRHG.measure_hom("p")
#     outcomes = CVRHG.hom_outcomes()
#     plt.figure(figsize=(16, 9))
#     plt.hist(outcomes, bins=100)
#     CVRHG.draw(label="hom_val_p")
#     CVRHG.draw(label="hom_val_p")
