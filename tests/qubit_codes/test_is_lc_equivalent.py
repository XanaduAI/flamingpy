"""Unit tests for is_lc_equivalent() method of EGraph class."""

import pytest
import numpy as np
from flamingpy.codes.graphs import EGraph
from flamingpy.utils import graph_states


def get_adj_mat(graph):
    """A convenience function to get the adjacency matrix of an EGraph object
    Args:
        graph: EGraph() object
    Returns:
        graph.adj_mat: numpy array corresponding to adjacency matrix
    """
    graph.adj_generator(sparse=False)
    return graph.adj_mat


def clifford_global_to_blocks(clifford_global):
    """A convenience function to get a dictionary of A, B, C, D blocks of
    clifford in global form:

    [A | B]
    [C | D]
    Args:
        clifford_global: a 2nx2n numpy array
    Returns:
        dict("A" : A, "B" : B, "C" : C, "D", D)
    """
    n = int(np.shape(clifford_global)[0] / 2)
    A = clifford_global[0:n, 0:n]
    B = clifford_global[0:n, n : 2 * n]
    C = clifford_global[n : 2 * n, 0:n]
    D = clifford_global[n : 2 * n, n : 2 * n]
    return {"A": A, "B": B, "C": C, "D": D}


def clifford_tensors_to_blocks(clifford_tensors):
    """A convenience function to get a dictionary of A, B, C, D blocks from
    clifford tensor factors:

    [A | B]
    [C | D]
    Args:
        clifford_tensors: a list of 2x2 numpy arrays
    Returns:
        dict("A" : A, "B" : B, "C" : C, "D", D)
    """
    n = int(len(clifford_tensors))
    A = np.diag([clifford_tensors[i][0, 0] for i in range(n)])
    B = np.diag([clifford_tensors[i][0, 1] for i in range(n)])
    C = np.diag([clifford_tensors[i][1, 0] for i in range(n)])
    D = np.diag([clifford_tensors[i][1, 1] for i in range(n)])
    return {"A": A, "B": B, "C": C, "D": D}


def mod2(M):
    """A convenience function that converts array elements to 0, 1 (modulo 2)
    Args:
        M : numpy array of any shape
    Returns:
        numpy array with elements in 0, 1
    """
    return np.vectorize(lambda x: x % 2)(M)


# computes matrix product and sum (modulo 2): X*Y+Z
def XYplusZ(X, Y, Z):
    """A convenience function to compute matrix product and sum (modulo 2): X*Y+Z
    Args:
        X: numpy array
        Y: numpy array
        Z: numpy array
    Returns:
        X*Y+Z: numpy array
    """
    return mod2(np.add(np.dot(X, Y), Z))


def distinct_tuples(min_range, max_range):
    """A convenience function that generates distinct tuples of integers in
    specified range.

    For example, can include (1,2) and (2,1) but not (2,2).
    Args:
        min_range (int): minimum of range
        max_range (int): maximum of range
    Returns:
        (list): list of (int, int) tuples
    """
    return [
        (a, b) for a in range(min_range, max_range) for b in range(min_range, max_range) if a != b
    ]


# The following tests are for the method 'is_lc_equivalent()' of the EGraph class.
# This method checks if two graphs are LC equivalent.
#
# The call signature of this method is:
# graph1.is_lc_equivalent(graph2, clifford_form) --> (equiv, clifford)
# Here graph1 and graph2 are EGraph objects,
# clifford_form is a string (either 'global' or 'tensor') which specifies the output form of the
# clifford if the equivalence is True.
#
# The method returns a tuple of the form (equiv, clifford).
# If the graphs are equivalent, 'equiv' is True and 'clifford' is the corresponding local clifford.
# If the graphs are not equivalent, 'equiv' is False and 'clifford' is None.
#
# Note the following asymmetry if graph1 and graph2 are LC equivalent:
# We can compare the graphs in two ways:
# graph1 --> graph2 :: by calling graph1.is_lc_equivalent(graph2)
# graph2 --> graph1 :: by calling graph2.is_lc_equivalent(graph1)
# Both return True equivalence, but the resulting cliffords may be distinct.
# Therefore, below we define test functions for both comparisons in order to test both cases.
#
# TESTING METHODOLOGY:
# When two graphs are equivalent, in order to verify the returned local clifford we perform the
# following algebraic check:
# Let G1 and G2 be the adjacency matrices of the two LC equivalent graphs, and let the local
# clifford operation be given in block form as:
# [A | B]
# [C | D]
# Then the following matrix block equation must hold True:
# G2*(C*G1+D) == (A*G1+B)
#
# Parametrize all tests to run for each of the two modes: "global" and "tensor"
@pytest.mark.parametrize("mode", ["global", "tensor"])
class TestLCEquivalent:
    """Test class for is_lc_equivalent() method of EGraph class."""

    # SKIPPED TEST
    # Runs 2 tests
    @pytest.mark.skip(
        reason="Current implementation assumes that equivalence with empty graph is False"
        "instead of raising a ValueError, so this test is irrelevant."
    )
    def test_emptygraph_with_emptygraph_assume_valueerror(self, mode):
        """Test if emptygraph equivalence raises a ValueError."""
        # Assert:
        with pytest.raises(ValueError):
            # ARRANGE:
            # define empty graph
            emptygraph = EGraph()
            # ACT:
            emptygraph.is_lc_equivalent(emptygraph, clifford_form=mode)

    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 2 tests
    def test_emptygraph_with_emptygraph_assume_notequivalent(self, mode):
        """Test if emptygraph equivalence returns (False, None)"""
        # ARRANGE:
        # define empty graph
        emptygraph = EGraph()
        # ACT:
        equiv, clifford = emptygraph.is_lc_equivalent(emptygraph, clifford_form=mode)
        # ASSERT:
        assert equiv is False
        assert clifford is None

    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 2 tests
    def test_pathgraph3_with_completegraph3_equivalent(self, mode):
        """Test True equivalence between path graph --> complete graph defined
        on 3 nodes."""
        # ARRANGE:
        # define all possible edges
        edge_1 = {(1, 0, 0), (0, 1, 0)}
        edge_2 = {(0, 1, 0), (0, 0, 1)}
        edge_3 = {(1, 0, 0), (0, 0, 1)}
        # construct path graph
        pathgraph3 = EGraph()
        pathgraph3.add_edges_from([edge_1, edge_2])
        # construct complete graph
        completegraph3 = EGraph()
        completegraph3.add_edges_from([edge_1, edge_2, edge_3])
        # get adjacency matrices
        adj1 = get_adj_mat(pathgraph3)
        adj2 = get_adj_mat(completegraph3)
        # ACT:
        # check equivalence pathgraph3 --> completegraph3
        equiv, clifford = pathgraph3.is_lc_equivalent(completegraph3, clifford_form=mode)
        # ASSERT:
        # get blocks of clifford depending on the form of output clifford
        if mode == "global":
            blocks = clifford_global_to_blocks(clifford)
        if mode == "tensor":
            blocks = clifford_tensors_to_blocks(clifford)
        # compute right-hand-side: (A*adj1 + B)
        rhs = XYplusZ(blocks["A"], adj1, blocks["B"])
        # compute left-hand-side: adj2*(C*adj1 + D)
        lhs = mod2(np.dot(adj2, XYplusZ(blocks["C"], adj1, blocks["D"])))
        assert equiv is True
        assert np.array_equal(lhs, rhs)

    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 2 tests
    def test_completegraph3_with_pathgraph3_equivalent(self, mode):
        """Test True equivalence between complete graph --> path graph defined
        on 3 nodes."""
        # ARRANGE:
        # define all possible edges
        edge_1 = {(1, 0, 0), (0, 1, 0)}
        edge_2 = {(0, 1, 0), (0, 0, 1)}
        edge_3 = {(1, 0, 0), (0, 0, 1)}
        # construct path graph
        pathgraph3 = EGraph()
        pathgraph3.add_edges_from([edge_1, edge_2])
        # construct complete graph
        completegraph3 = EGraph()
        completegraph3.add_edges_from([edge_1, edge_2, edge_3])
        # get adjacency matrices
        adj1 = get_adj_mat(completegraph3)
        adj2 = get_adj_mat(pathgraph3)
        # ACT:
        # check equivalence completegraph3 --> pathgraph3
        equiv, clifford = completegraph3.is_lc_equivalent(pathgraph3, clifford_form=mode)
        # ASSERT:
        # get blocks of clifford depending on the form of output clifford
        if mode == "global":
            blocks = clifford_global_to_blocks(clifford)
        if mode == "tensor":
            blocks = clifford_tensors_to_blocks(clifford)
        # compute right-hand-side: (A*adj1 + B)
        rhs = XYplusZ(blocks["A"], adj1, blocks["B"])
        # compute left-hand-side: adj2*(C*adj1 + D)
        lhs = mod2(np.dot(adj2, XYplusZ(blocks["C"], adj1, blocks["D"])))
        assert equiv is True
        assert np.array_equal(lhs, rhs)

    # Tests graphs on same number of nodes in range(1,5)
    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 8 tests
    @pytest.mark.parametrize("nodes", range(1, 5))
    def test_completegraph_with_stargraph_equivalent(self, nodes, mode):
        """Test True equivalence between complete graph --> star graph defined
        on same number of nodes."""
        # ARRANGE:
        # define complete graph
        graph1 = graph_states.complete_graph(nodes)
        # define star graph
        graph2 = graph_states.star_graph(nodes)
        # get adjacency matrices
        adj1 = get_adj_mat(graph1)
        adj2 = get_adj_mat(graph2)
        # ACT:
        # check equivalence complete graph --> star graph
        equiv, clifford = graph1.is_lc_equivalent(graph2, clifford_form=mode)
        # ASSERT:
        # get blocks of clifford depending on the form of output clifford
        if mode == "global":
            blocks = clifford_global_to_blocks(clifford)
        if mode == "tensor":
            blocks = clifford_tensors_to_blocks(clifford)
        # compute right-hand-side: (A*adj1 + B)
        rhs = XYplusZ(blocks["A"], adj1, blocks["B"])
        # compute left-hand-side: adj2*(C*adj1 + D)
        lhs = mod2(np.dot(adj2, XYplusZ(blocks["C"], adj1, blocks["D"])))
        assert equiv is True
        assert np.array_equal(lhs, rhs)

    # Tests graphs on same number of nodes in range(1,5)
    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 8 tests
    @pytest.mark.parametrize("nodes", range(1, 5))
    def test_stargraph_with_completegraph_equivalent(self, nodes, mode):
        """Test True equivalence between star graph --> complete graph defined
        on same number of nodes."""
        # ARRANGE:
        # interchange graph1 and graph2 compared to dual test
        # define complete graph
        graph2 = graph_states.complete_graph(nodes)
        # define star graph
        graph1 = graph_states.star_graph(nodes)
        # get adjacency matrices
        adj1 = get_adj_mat(graph1)
        adj2 = get_adj_mat(graph2)
        # ACT:
        # check equivalence star graph --> complete graph
        equiv, clifford = graph1.is_lc_equivalent(graph2, clifford_form=mode)
        # ASSERT:
        # get blocks of clifford depending on the form of output clifford
        if mode == "global":
            blocks = clifford_global_to_blocks(clifford)
        if mode == "tensor":
            blocks = clifford_tensors_to_blocks(clifford)
        # compute right-hand-side: (A*adj1 + B)
        rhs = XYplusZ(blocks["A"], adj1, blocks["B"])
        # compute left-hand-side: adj2*(C*adj1 + D)
        lhs = mod2(np.dot(adj2, XYplusZ(blocks["C"], adj1, blocks["D"])))
        assert equiv is True
        assert np.array_equal(lhs, rhs)

    # Tests graphs on different number of nodes in range(1,5)
    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 24 tests
    @pytest.mark.parametrize("nodes1, nodes2", distinct_tuples(1, 5))
    def test_completeraph_with_stargraph_notequivalent(self, nodes1, nodes2, mode):
        """Test False equivalence between complete graph --> star graph defined
        on different number of nodes."""
        # ARRANGE:
        # define complete graph
        graph1 = graph_states.complete_graph(nodes1)
        # define star graph
        graph2 = graph_states.star_graph(nodes2)
        # ACT:
        # check equivalence complete graph --> star graph
        equiv, clifford = graph1.is_lc_equivalent(graph2, clifford_form=mode)
        # ASSERT:
        assert equiv is False
        assert clifford is None

    # Tests graphs on different number of nodes in range(1,5)
    # Tests both output clifford_form modes: 'global' and 'tensor'
    # Runs 24 tests
    @pytest.mark.parametrize("nodes1, nodes2", distinct_tuples(1, 5))
    def test_stargraph_with_completegraph_notequivalent(self, nodes1, nodes2, mode):
        """Test False equivalence between star graph --> complete graph defined
        on different number of nodes."""
        # ARRANGE:
        # interchange graph1 and graph2 compared to dual test
        # define complete graph
        graph2 = graph_states.complete_graph(nodes1)
        # define star graph
        graph1 = graph_states.star_graph(nodes2)
        # ACT:
        # check equivalence star graph --> complete graph
        equiv, clifford = graph1.is_lc_equivalent(graph2, clifford_form=mode)
        # ASSERT:
        assert equiv is False
        assert clifford is None
