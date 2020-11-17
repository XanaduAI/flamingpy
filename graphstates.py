# Copyright 2020 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
'''This module provides classes for DV and CV graph states.'''
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Style for figures
fontdict = {'fontsize': 14, 'family': 'serif'}


class EGraph(nx.Graph):
    '''An enhanced graph class based on networkx.Graph.'''

    def __init__(self):
        nx.Graph.__init__(self)
        if 'dims' not in self.graph:
            self.graph['dims'] = None

    def adj(self):
        return nx.to_numpy_array(self)

    def color(self, coord):
        return self.nodes[coord]['color']

    def draw(self, color=1, label=0):
        '''Draw the graph.
        Args:
            label (bool): if True, label the indices; unlabelled by
                default.
            color (bool): if True, color the syndrome and data qubits
                red and green, respectively; colored by default.
        Returns:
            None
        '''
        # Recommended to be viewed with IPython.
        if self.graph['dims']:
            nx, ny, nz = self.graph['dims']
        else:
            nx, ny, nz = 5, 5, 5
        fig = plt.figure(figsize=(2 * (nx + ny + nz + 2),
                                  2 * (nx + ny + nz + 2)))
        ax = fig.add_subplot(111, projection='3d')
        # Plotting points. y and z are swapped in the loops to compare
        # the lattice diagrams from our notes, where the z axis goes
        # into the page; however the axes labels are correct.
        for point in self.nodes:
            x, z, y = point
            node_color = self.color(point)
            ax.scatter(x, y, z, s=70, c=color*node_color+(1-color)*'k')
            indices = {c: n for (n, c) in enumerate(self.nodes)}
            if label:
                ax.text(x, y, z, str(indices[point]), fontdict=fontdict,
                        color='MediumBlue', backgroundcolor='w')
        # Plotting edges.
        for edge in self.edges:
            x1, z1, y1 = edge[0]
            x2, z2, y2 = edge[1]
            plt.plot([x1, x2], [y1, y2], [z1, z2], c='grey')

        plt.xticks(range(0, 2*nx + 1))
        plt.yticks(range(0, 2*nz + 1))
        ax.set_zticks(range(0, 2*ny + 1))
        ax.set_xlabel('x', fontdict=fontdict)
        ax.set_ylabel('z', fontdict=fontdict)
        ax.set_zlabel('y', fontdict=fontdict)
        plt.rcParams['grid.color'] = "lightgray"
        plt.tight_layout(pad=3)
        plt.draw()
        return ax

    def index(self):
        indexed_graph = nx.convert_node_labels_to_integers(self, ordering='sorted', label_attribute='pos')
        return indexed_graph


def RHG_graph(dims, pol=0):
    '''Return an EGraph of the RHG lattice.'''
    # Dimensions of the lattice.
    if np.size(dims) == 1:
        dims = (dims, dims, dims)
    nx, ny, nz = dims

    lattice = EGraph()
    lattice.graph['dims'] = dims
    # Coordinates of red qubits in even and odd vertical slices.
    even_red = [(2*i+1, 2*j+1, 2*k) for (i, j, k) in
                it.product(range(nx), range(ny), range(nz+1))]
    odd_red = [(2*i, 2*j, 2*k+1) for (i, j, k) in
               it.product(range(nx+1), range(ny+1), range(nz))]
    all_red = set(even_red + odd_red)

    # Coordinates of green qubits in even and odd horizontal slices.
    even_green = [(2*i+1, 2*j, k) for (i, j, k) in
                  it.product(range(nx), range(ny+1), range(2*nz+1))]
    odd_green = [(2*i, 2*j+1, k) for (i, j, k) in
                 it.product(range(nx+1), range(ny), range(2*nz+1))]
    all_green = set(even_green + odd_green)

    # Coordinates of all potential neighbours of red vertices.
    def red_neighbours(p):
        right = (p[0]+1, p[1], p[2])
        top = (p[0], p[1]+1, p[2])
        left = (p[0]-1, p[1], p[2])
        bottom = (p[0], p[1]-1, p[2])
        return [bottom, left, top, right]

    def green_neighbours(p):
        return {(p[0], p[1], p[2]-1), (p[0], p[1], p[2]+1)}

    for point in all_red:
        for i in range(4):
            polarity = (-1) ** (pol * (point[2] + i))
            neighbour = red_neighbours(point)[i]
            if neighbour in all_green:
                lattice.add_edge(point, neighbour, weight=polarity)
        lattice.nodes[point]['color'] = 'red'

    for point in all_green:
        polarity = (-1) ** (pol * (point[1] + 1))
        for neighbour in green_neighbours(point):
            if neighbour in all_green:
                lattice.add_edge(point, neighbour, weight=polarity)
        lattice.nodes[point]['color'] = 'green'

    return lattice


if __name__ == '__main__':
    R = RHG_graph(1)
    R.draw(label=1)
    R1 = R.index()
    A = R.adj()
    A1 = R1.adj()
