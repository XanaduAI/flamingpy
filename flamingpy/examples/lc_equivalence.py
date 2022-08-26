"""Example for testing LC equivalence of graph states."""
from flamingpy.utils.graph_states import star_graph, complete_graph, linear_cluster, ring_graph

print("Testing LC equivalence of graph states:", "\n")

# number of graph nodes (qubits)
nodes = 3
print("Number of nodes/qubits: ", nodes)

# output form of local Clifford (can be 'tensor' or 'global')
clifford_form = "tensor"  # 'global
print("Clifford form: ", f"{clifford_form}", "\n")

# test LC equivalence with star_graph and complete_graph
# equivalent for any number of nodes
equiv, clifford = star_graph(nodes).is_lc_equivalent(complete_graph(nodes), clifford_form)
print("star_graph --> complete_graph:")
print("Equivalent: ", equiv)
print(f"Local Clifford ({clifford_form}): ", clifford, "\n")

# test LC equivalence with star_graph and linear_cluster
# equivalent only if number of nodes <= 3
equiv, clifford = star_graph(nodes).is_lc_equivalent(linear_cluster(nodes), clifford_form)
print("star_graph --> linear_cluster")
print("Equivalent: ", equiv)
print(f"Local Clifford ({clifford_form}): ", clifford, "\n")

# test LC equivalence with star_graph and ring_graph
# equivalent only if number of nodes <= 3
equiv, clifford = star_graph(nodes).is_lc_equivalent(ring_graph(nodes), clifford_form)
print("star_graph --> ring_graph:")
print("Equivalent: ", equiv)
print(f"Local Clifford ({clifford_form}): ", clifford, "\n")
