import networkx as nx
import numpy as np
import time
from ft_stack.lemon import max_weight_matching
import matplotlib.pyplot as plt


time_lemon = []
time_nx = []
matrix_size = []

# Change that for more data points
for num_nodes in range(50,150,50):

    matrix_size.append(num_nodes)
    b = np.random.random_integers(0,2000,size=(num_nodes, num_nodes))
    b_symm = (b + b.T)*.5
    for i in range(num_nodes):
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
plt.savefig("lemon_benchmark.pdf")

