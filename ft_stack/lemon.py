import lemonpy
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import time
import sys

def max_weight_matching(G_match, weight):
    """Compute the maximum weighted matching graph using lemon.
       Assumptions:
          1. Symmetric adjacency matrix.
          2. Adjacency matrix has zeros along diagonal.
    """
    adjacency = nx.to_numpy_matrix(G_match,weight=weight)
    lemon_matching = lemonpy.mwpm(adjacency)

    #lemon uses different node ids, so we convert back to networkx node ids
    nx_map = list(G_match.nodes())
    matching = set()
    for i in lemon_matching:
        matching.add((nx_map[i[0]],nx_map[i[1]]))
    return matching


if __name__ == "__main__":
    
    time_lemon = []
    time_nx = []
    matrix_size = []
    
    for N in range(50,550,50):

        matrix_size.append(N)
        b = np.random.random_integers(0,2000,size=(N,N))
        b_symm = (b + b.T)*.5
        for i in range(N):
            b_symm[i,i] = 0
        
        G_match=nx.from_numpy_matrix(b_symm)
        
        start = time.time()
        lemon_matching = max_weight_matching(G_match,weight='weight')
        end = time.time()
        time_lemon.append(end-start)

        start = time.time()
        nx_matching = nx.max_weight_matching(G_match,weight='weight')
        end = time.time()
        time_nx.append(end-start)

        print("matrix_size t_lemon t_nx = ",matrix_size[-1], time_lemon[-1], time_nx[-1])
    
    plt.plot(matrix_size, time_lemon,label="lemon")
    plt.plot(matrix_size, time_nx, label="networkx")
    plt.xlabel("Nodes")
    plt.ylabel("Time")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig("benchmark.png")
    
