import ft_stack.lemonpy as lp
import networkx as nx


def max_weight_matching(G_match, weight):
    """Compute the maximum weighted matching graph using lemon.
       Assumptions:
          1. Symmetric adjacency matrix.
          2. Adjacency matrix has zeros along diagonal.
    """
    adjacency = nx.to_numpy_matrix(G_match,weight=weight)
    lemon_matching = lp.mwpm(adjacency)

    #lemon uses different node ids, so we convert back to networkx node ids
    nx_map = list(G_match.nodes())
    matching = set()
    for i in lemon_matching:
        matching.add((nx_map[i[0]],nx_map[i[1]]))
    return matching
