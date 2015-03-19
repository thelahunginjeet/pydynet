#!/usr/bin/env python
# encoding: utf-8
"""
rewiring.py

Methods for growing/shrinking/rewiring graphs.  These functions modify the input graphs in-place!

Created by Kevin Brown on 2014-12-04.
"""

import networkx as nx
from utilities import randchoice


def add_random_edge(G):
    """
    Adds an edge between two nodes.  The first node is randomly chosen, and then another
    node is randomly chosen from the set of nodes not already connected to the first node.
    Nothing happens if the first node selected is connected to every other node in the graph.

    INPUT:
    -------
    G : networkx graph, required
        changed in-place!

    OUTPUT:
    -------
    flag : boolean
        indicates if the add was successful or not

    -Conservative in nodes, nonconservative in edges.
    -Changes degree distribution.
    -Fails if the first node selected is fully connected.
    """
    # pick a random node
    n1 = randchoice(G.nodes())
    # neighbors of n1
    nb = G.neighbors(n1)
    validNodes = [x for x in G.nodes() if x != n1 and x not in nb]
    if len(validNodes) == 0:
        return False
    # second unique node, not already connected to the first node
    n2 = randchoice(validNodes)
    G.add_edge(n1,n2)
    return True


def remove_random_edge(G):
    """
    Removes a randomly chosen edge in the graph.

    INPUT:
    -------
    G : networkx graph, required
        changed in-place!

    OUTPUT:
    -------
    flag : boolean
        indicates if the delete was successful or not

    -Conservative in nodes, nonconservative in edges.
    -Changes degree distribution.
    -Fails for a graph with zero edges.
    """
    # choose existing edge at random and unpack tuple
    validEdges = G.edges()
    if len(validEdges) == 0:
        print 'Graph has no edges to remove!'
        return False
    e = randchoice(validEdges)
    G.remove_edge(*e)
    return True


def move_random_edge(G):
    """
    Chooses a random edge and a random additional node, not contained in that edge, and
    moves the edge to terminate at the new node.  This operation might disconnect the
    graph (leave singleton nodes).

    INPUT:
    -------
    G : networkx graph, required
        changed in-place!

    OUTPUT:
    -------
    flag : boolean
        indicates if the move was successful or not

    -Conservative in nodes and edges.
    -Changes degree distribution.
    -Fails if the first node selected is fully connected or for a graph with no edges.
    """
    # choose a random edge
    validEdges = G.edges()
    if len(validEdges) == 0:
        print 'Graph has no edges to move!'
        return False
    e = randchoice(validEdges)
    # neighbors of the chosen edge
    nb = G.neighbors(e[0])
    # additional node that is neither n1 or connected to n1
    validNodes = [x for x in G.nodes() if x != e[0] and x not in nb]
    if len(validNodes) == 0:
        return False
    n3 = randchoice(validNodes)
    G.remove_edge(*e)
    G.add_edge(e[0],n3)
    return True


def swap_random_edges(G):
    """
    Exchanges termnii of two edges e1 and e2, such that:
        (n1,n2); (n3,n4) --> (n1,n4); (n3,n2)

    INPUT:
    -------
    G : networkx graph, required
        changed in-place!

    OUTPUT:
    -------
    flag : boolean
        indicates if the swap was successful or not

    -Conservative in nodes and edges.
    -Does not change degree distribution. (Does not even change individual node degrees.)
    """
    # choose an edge at random
    validEdges = G.edges()
    if len(validEdges) < 2:
        print 'Graph does not have enough edges to swap!'
        return False
    n1,n2 = randchoice(validEdges)
    # choose a third node such that
    #   -n3 is not either n1 or n2
    #   -n3 is not a neighbor of n2
    nb = G.neighbors(n2)
    validNodes = [x for x in G.nodes() if x != n1 and x != n2 and x not in nb]
    if len(validNodes) == 0:
        return False
    n3 = randchoice(validNodes)
    # choose a fourth node such that
    #   -n4 is not either n1 or n2
    #   -n4 is connected to n3
    #   -n4 is not a neighbor of n1
    nb3 = G.neighbors(n3)
    nb1 = G.neighbors(n1)
    validNodes = [x for x in nb3 if x != n1 and x != n2 and x not in nb1]
    if len(validNodes) == 0:
        return False
    n4 = randchoice(validNodes)
    # do the swap
    G.remove_edges_from([(n1,n2),(n3,n4)])
    G.add_edges_from([(n1,n4),(n2,n3)])
    return True


