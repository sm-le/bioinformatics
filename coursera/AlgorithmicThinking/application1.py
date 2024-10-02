"""Algorithmic Thinking Application 1:
by smlee

Analysis of Citation Graphs
"""

"""
Question 1: Compute the in-degree distribution for given citation graph. The distribution should be normalized; and then, compute a log plot of the points.

Assessment criteria
- Does the plot follow the formatting guidelines for plots?
- Is the plot that of a normalized distribution on a log/log scale?
- Is the content of the plot correct?
"""

# Modules
import urllib
from typing import Dict, Set
import matplotlib.pyplot as plt

# Main
def get_graph_from_url(url:str | None=None) -> Dict[int,Set[int]]:
    """Request data from input url and construct a graph

    Args:
        url: a resource locator
    Returns:
        graph in a dictionary(int,set) format
    """
    # request data from url
    with urllib.request.urlopen(url) as conn:
        text = list()
        while True:
            line = conn.readline().decode().strip()
            if not line:
                break
            text.append(line)
    
    # construct graph
    graph = dict()
    for nodeEdge in text:
        components = list(map(lambda x: int(x), nodeEdge.split(" ")))
        node, edges = components[0], set(components[1:]) if len(components) > 1 else set()
        if node not in graph:
            graph[node] = set()
        graph[node].update(edges)

    return graph

def compute_in_degrees(digraph:dict | None=None) -> Dict[int,int]:
    """Compute in-degrees for the nodes given digraph (dictionary)

    Args:
        digraph: a dictionary representation of a graph
    Returns:
        in-degrees graph in a dictionary
    """
    if digraph is None:
        raise ValueError("Input is missing. Please use digraph")
    # get empty graph
    graph = {i:0 for i in digraph.keys()}

    # get all directional edges from digraph and loop to count
    for node in digraph.values():
        for edge in node:
            graph[edge] += 1
    
    return graph

def normalized_in_degree_distribution(digraph:dict) -> Dict[int,int]:
    """Compute normalized in-degrees distribution given digraph

    Args:
        digraph: a dictionary representation of a graph
    Returns:
        Normalized in-degrees distribution in a dictionary
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

    # calculate sum
    total = sum(distbn.values())
    # normalize
    distbn = dict(map(lambda x: (x[0],x[1]/total), distbn.items()))
    
    return distbn

def plot_graph(data:dict,
               title:str,
               *,
               log:bool=True,
               fn:str=None) -> plt.scatter:
    """Plot scatterplot using x data

    Args:
        data: data for plotting
        log: do log scale transformation?
        fn: filename
    Returns:
        plt.scatter plot
    """

    x = data.keys()
    y = data.values()
    
    fig, ax = plt.subplots(figsize=(9,6))
    
    ax.scatter(x,y, s=60, alpha=0.7, edgecolors="k")
    ax.set_xlabel("Number of citations")
    ax.set_ylabel("Normalized in-degree distribution")
    ax.set_title(title)
    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")

    if fn:
        plt.savefig(fn)




"""
Question 2: Adjust algorithm of make_complete_graph function to add edge with probability p

Assessment criteria:
- Is the expected value of the in-degree the same for every node in an ER graph? Please answer yes or no and include a short explanation for your answer.
    -> yes, if there are n nodes and probability p, the expected in-degree of each node is (n-1)p since they are equally likely to receive an edge from other nodes.
- What does the in-degree distribution for an ER graph look like? Provide a short written description of the shape of the distribution.
    -> It is binomially distibuted with centered around the expected value.
- Does the shape of the in-degree distribution plot for ER look similar to the shape of the in-degree distribution for the citation graph? Provide a short explanation of the similarities or differences. Focus on comparing the shape of the two plots as discussed in the class page on "Creating, formatting, and comparing plots".
    -> The shape of the in-degree distribution plot for ER looks different to the shape of the in-degree distribution for the citation graph. The shape of the citation graph is slightly linear with falling trend to the right while the ER graph is randomly distributed with rising trend toward the expected value and falls right after. 
"""

import random
def make_stochastic_graph(num_nodes:int | None=None,
                          *,
                          p:float=0.5) -> Dict[int,Set[int]]:
    """Make a stochastic graph given a number of nodes

    Args:
        num_nodes (int): number of nodes
        p (float) = probability, default: 0.5
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
        graph[count] = {i for i in range(num_nodes) if random.random() > p} - {count}
        count += 1

    return graph

"""
Question 3 & 4: Synthetic graph generation with n and m

Assessment criteria:
    Provide numerical values for n and m
    Provide plot for the synthetic graphfirefox
"""

def make_complete_graph(num_nodes:int | None=None) -> Dict[int,Set[int]]:
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

def make_DPA_graph(num_nodes:int,
                   m:int) -> Dict[int,Set[int]]:
    """Make a probabilistic graph using DPA algorithm

    Args:
        num_nodes (int): number of nodes
        m (int): number of fixed number of nodes to start with
    Returns:
        directed graph in a dictionary
    """

    # initialize a graph with m nodes
    graph = make_complete_graph(m)
    # initial graph summary
    in_degree = compute_in_degrees(graph)
    current_nodes = list(in_degree.keys())
    current_node_size = len(current_nodes)
    total_in_degree = sum(in_degree.values())
    # loop
    for i in range(m,num_nodes): 
        # select nodes with probability p
        p = [(in_degree[j] + 1) / (total_in_degree + current_node_size) for j in range(current_node_size)]
        edges = set(random.choices(current_nodes, weights=p, k=m))
        # update graph
        graph[i] = edges
        current_nodes.append(i)
        current_node_size += 1
        total_in_degree += len(edges)
        # update in-degree
        for e in edges:
            in_degree[e] += 1
            in_degree[i] = 0
    return graph
        
"""
Question 5: Comparison between DPA and citation graph

The in-degree distribution for the DPA graph is similar to the citation graph. We assigned drawing probability such that nodes with high 
number of edges are more likely to be sampled during the process. This phenomena is close representation of rich gets richer phenomenon as 
highly connected nodes are more visible; and hence, gets more citation. Even in physics citation papers, researchers will favor highly cited 
paper unless they have bias. Therefore, rich get richers phenomenon applies here. 
"""

if __name__ == "__main__":
    # Constant
    URL_CITATION = "https://storage.googleapis.com/codeskulptor-alg/alg_phys-cite.txt"

    # Question 1
    g = get_graph_from_url(URL_CITATION)
    dist = normalized_in_degree_distribution(g)
    plot_graph(dist,title="in-degree distribution of the citation graph",fn="inDegreeDistribution.png")

    # Question 2
    NO_NODES = 10
    s = make_stochastic_graph(NO_NODES)
    sdist = normalized_in_degree_distribution(s)
    plot_graph(sdist) 

    # Question 3 & 4
    g = make_DPA_graph(27770,13)
    ng = normalized_in_degree_distribution(g)
    plot_graph(ng,title="in-degree distribution of the DPA graph",fn="inDegreeDistributionDpaGraph.png")