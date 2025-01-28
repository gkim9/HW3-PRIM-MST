import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import numpy as np
from mst.graph import Graph
from sklearn.metrics import pairwise_distances

def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001,
              ):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    total = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
            total += mst[i, j]
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'

    # number of edges in a minimum spanning tree is always equal to total number of nodes - 1 (since cannot be cyclic and goes through all nodes)
    # counting the instances of all non-zero edges in mst and dividing by 2 since edges are double counted (i -> j == j -> i) to get number of edges in mst
    assert np.count_nonzero(mst)/2 == mst.shape[0] - 1, "Incorrect amount of edges in minimum spanning tree."

def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)

print(test_mst_small())

def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = './data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)

def test_mst_student():
    """
    
    TODO: Write at least one unit test for MST construction.

    Testing if the input graph has an edge from a node to iself with a weight not equal to 0.
    
    """

    # Input graph has a value of non-zero in the diagonal
    graph = Graph('data/non_zero_diagonal.csv')
    with pytest.raises(Exception):
        graph.construct_mst()

    # Input graph is not connected -- should not have the right number of edges to be a minimum spanning tree
    not_connected = Graph('data/non_connected_example.csv')
    not_connected.construct_mst()

    with pytest.raises(AssertionError):
        check_mst(not_connected.adj_mat, not_connected.mst, 9)
