import numpy as np
from flamingpy.codes.graphs import EGraph
from flamingpy.utils import graph_states

graph1 = graph_states.complete_graph(4)

graph2 = graph_states.star_graph(4)

output = graph1.is_lc_equivalent(graph2, 'tensor')

print(output)