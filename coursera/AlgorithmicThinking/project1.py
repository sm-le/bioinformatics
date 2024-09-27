"""Algorithmic Thinking Project 1:
by smlee

Python Functions for Degree Distribution for Graphs
"""

# Constants
EX_GRAPH0 = {0:set([1,2]),
             1:set(),
             2:set()}
EX_GRAPH1 = {0:set([1,4,5]),
             1:set([2,6]),
             2:set([3]),
             3:set([0]),
             4:set([1]),
             5:set([2]),
             6:set()}
EX_GRAPH2 = {0:set([1,4,5]),
             1:set([2,6]),
             2:set([3,7]),
             3:set([7]),
             4:set([1]),
             5:set([2]),
             6:set(),
             7:set([3]),
             8:set([1,2]),
             9:set([0,3,4,5,6,7])}

# Modules
# from typing import Dict, Set
# Python functions
def make_complete_graph(num_nodes):#:int):# | None=None) -> Dict[int,Set[int]]:
    """Make a complete graph given a number of nodes

    Args:
        num_nodes (int): number of nodes
    Returns:
        directed graph in a dictionary
    """
    # set initial value
    graph = dict()
    # check input
    if num_nodes is None:
        return graph 

    # set initial node number
    count = 0
    # while to iterate and add node to graph
    while num_nodes > 0 \
        and count < num_nodes:
        graph[count] = {i for i in range(num_nodes)} - {count}
        count += 1

    return graph

def compute_in_degrees(digraph):#:dict):# -> Dict[int,Set[int]]:
    """Compute in-degrees for the nodes given digraph (dictionary)

    Args:
        digraph: a dictionary representation of a graph
    Returns:
        in-degrees graph in a dictionary
    """
    # get empty graph
    graph = {i:0 for i in digraph.keys()}

    # get all directional edges from digraph and loop to count
    for node in digraph.values():
        for edge in node:
            graph[edge] += 1
    
    return graph

def in_degree_distribution(digraph):#:dict):# -> Dict[int,int]:
    """Compute in-degrees distribution for the nodes given digraph (dictionary)

    Args:
        digraph: a dictionary representation of a graph
    Returns:
        in-degrees graph in a dictionary

    """

    # get in-degrees graph
    graph = compute_in_degrees(digraph)

    # set empty distribution
    distbn = dict()
    # loop graph to count no. in-degrees 
    for in_degree in graph.values():
        if in_degree not in distbn:
            distbn[in_degree] = 0
        distbn[in_degree] += 1

    return distbn