![BuildStatus](https://github.com/gkim9/HW3-PRIM-MST/actions/workflows/test.yml/badge.svg?event=push)

# HW 3: Prim's algorithm

## Description
Prim's algorithm can find the Minimum Spanning Tree (MST) of a graph, provided that the graph is fully connected and undirected. MST is a graph such that all vertices are connected with no cycles.

The algorithm begins by picking an arbitrary starting vertex and by choosing the smallest edge from it, exploring the connected vertex. The newly explored vertex is added to the set of explored vertices, and the smallest edge from any of the vertices in the explored set to a vertex not in the explored set is added. This process continues until all nodes have been explored. 

mst/graph.py will create an "mst" attribute of a numpy array containing the weights of the edges of the mst. To achieve this, we use the heapq module to create a sorted queue.

## Provided README
In this assignment, you'll implement Prim's algorithm, a non-trivial greedy algorithm used to construct minimum spanning trees. 

## Tasks

### Coding

- [x] Complete the `construct_mst` method found in `mst/graph.py`. All necessary modules have already been imported. You should not rely on any other dependencies (e.g. networkx). 

### Development

- [x] Add more assertions to the `check_mst` function in `test/test_mst.py`.
- [x] Write at least one more unit test (in the `test_mst.py` file) for your `construct_mst` implementation. (Two unit tests have already been provided: the first operates on a small graph of four nodes, and the second on a larger graph of 140 single cells, projected onto a lower dimensional subspace.)
- [x] Make your package `pip` installable. (Refer to prevous assignments for more in-depth information.)
- [x] Automate testing with `pytest` and GitHub Actions, and add a status badge to this README file. (Refer to previous assignments for more in-depth information.)

## Getting started

Fork this repository to your own GitHub account. Work on the codebase locally and commit changes to your forked repository. 

You will need following packages:

- [numpy](https://numpy.org/)
- [scikit-learn](https://scikit-learn.org/)
- [pytest](https://docs.pytest.org/en/7.2.x/)

We also strongly recommend you use the built-in [heapq](https://docs.python.org/3/library/heapq.html) module.

## Grading

### Code (6 points)

* Minimum spanning tree construction works correctly (6)
    * Correct implementation of Prim's algorithm (4)
    * Produces expected output on small graph (1) 
    * Produces expected output on single cell data (1) 

### Unit tests (3 points)

* Complete function "check_mst" (1)
* Write at least two unit tests for MST construction (2)

### Style (1 points)

* Readable code with clear comments and method descriptions (1)

### Extra credit (0.5)

* Github actions/workflow (0.5)
