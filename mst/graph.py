import numpy as np
import heapq
from typing import Union

class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')
    
    def check_input(self):
        for node in range(self.adj_mat.shape[0]):
            print(range(self.adj_mat.shape[0]))
            if self.adj_mat[node][node] != 0:
                print(node, self.adj_mat[node][node])
                return "Error"

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.


        Pseudocode provided in class as a reference:
        initialize S and T, s is any node in set of vertices (V)
        For each v != s, pi[v] <- inf, pred[v] <- null, pi[s] <- 0      distance from s -> v is infinity, no predecessor to v, and distance to s from s is 0
        Create an empty priority queue (pq)

        For each vertex in V, insert v in priority queue (place sorted based on the lowest edge weight/cost between v and S)
        While pq is not empty,
            take first item in pq = u
            add u to the set S of "explored" nodes, add pred[u] into T
            for each edge e = (u, v) in Edge set with v not in the explored set (S),
                if cost of the edge e is less than the pi[v] (cheapest known edge between v and S)
                decrease-key of the priority queue of that vertext based on the new smallest edge
                the smallest edge becomes the new pi[v] and pred[v] becomes e

        """
        if self.check_input() == "Error":
            raise Exception(f"Tree is incorrect; weight along a diagonal is not 0")
        self.num_nodes = self.adj_mat.shape[0]
        node_list = list(range(self.num_nodes))
        # initializing the mst with all weights of edges = 0
        self.mst = np.zeros((self.num_nodes, self.num_nodes))

        explored_nodes = set()
        edges = []

        pi_dict = {node: float('inf') for node in range(self.num_nodes)}
        pred_dict = {node: None for node in range(self.num_nodes)}

        start_node = np.random.choice(list(range(self.num_nodes))) #choosing start node
        pi_dict[start_node] = 0 # reflecting pi[s] = 0
        pred_dict[start_node] = start_node # starting node's parent is itself

        priority_queue = []
        heapq.heapify(priority_queue) # creating priority queue that will be sorted
        for node in range(self.num_nodes):
            heapq.heappush(priority_queue, [pi_dict[node], node]) # placing the weight of the node in the front so that priority queue can be sorted by the weight
        
        while priority_queue:
            weight_u, node_u = heapq.heappop(priority_queue)

            # Skip if already explored
            if node_u in explored_nodes:
                continue 

            explored_nodes.add(node_u)

            if node_u != start_node and pred_dict[node_u] is not None:
                edges.append([pred_dict[node_u], node_u])

            for node_v, weight_v in enumerate(self.adj_mat[node_u]):
                if node_v not in explored_nodes and weight_v != 0:
                    if weight_v < pi_dict[node_v]:
                        pi_dict[node_v] = weight_v
                        pred_dict[node_v] = node_u
                        heapq.heappush(priority_queue, [pi_dict[node_v], node_v])

        for node_u, node_v in edges:
            self.mst[node_u][node_v] = self.adj_mat[node_u][node_v]
            self.mst[node_v][node_u] = self.adj_mat[node_u][node_v]
            print(edges)
        print("THISISIT", self.mst)

# example = Graph('data/small.csv')
# example.construct_mst()
# print(example.mst)
from sklearn.metrics import pairwise_distances
file_path = './data/slingshot_example.txt'
coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
g = Graph(dist_mat)
g.construct_mst()

tot = np.sum(g.mst)

print(g.mst)
print(tot)