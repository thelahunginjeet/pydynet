"""
network.py

This module controls the construction of a network of oscillators.

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

"""

from __future__ import division
from utilities import randchoice,randspin
from numpy.random import rand
from numpy import zeros,ones,arange,empty_like,where,reshape,asarray
from numpy import sqrt,cos,sin,pi,mod,round,all
from numpy import mean,var
from numpy import uint8,float64,int
import networkx as nx
import eulerint


def dydtMS(y,t,p):
    """
    Mirollo-Strogatz linear governing equation for oscillator dynamics.

            dydt = -p[0]*y + p[1]

    Parameter values from Mirollo and Strogatz 1989:
            p[0] = 1
            p[1] = 2
    """
    return -p[0]*y + p[1]



class PulseOscillatorNetwork(nx.Graph):

    def __init__(self,N,topology,*args):
        super(PulseOscillatorNetwork,self).__init__()
        # dispatch on topology type
        tdict = {'empty':self.connect_empty, 'full':self.connect_full, 'ring':self.connect_ring, 'fixed degree':self.connect_fixed_degree,
                 'fixed edges':self.connect_fixed_edges,'ERnp':self.connect_gnp, 'WS':self.connect_watts_strogatz,
                 'NWS':self.connect_newman_watts_strogatz,'BA':self.connect_barabasi_albert,'ERnm':self.connect_gnm,
                 'configuration':self.connect_configuration, 'edgelist':self.connect_edgelist}
        if tdict.has_key(topology):
            tdict[topology](N,*args)
        else:
            print('ERROR.  Unrecognized graph topology. Defaulting to ring.')
            self.connect_ring(N)
        # set default amplitude/delay/threshold parameters for dynamical simulations
        self.set_synaptic('excitatory')
        self.delta = 0
        self.y_th = 1.0
        # default distance embedding
        self.embed = DistanceEmbedding(N)
        self.embed.unitcirc_map()
        self.set_edge_lengths(self.embed)

    def set_synaptic(self,mode='excitatory',p=0.5):
        """
        Sets up the matrix of synaptic weights.  Currently supported options are:
            'excitatory': all connections have weight +1/(N-1)
            'inhibitory': all connections have weight -1/(N-1)
            'random' : connections are inhibitory with probability p and excitatory
                with probability 1-p
        """
        N = len(self.nodes())
        self.eps = (1.0/(N-1))*ones((N,N))
        if mode is 'inhibitory':
            self.eps = -1*self.eps
        elif mode is 'random':
            for e in self.edges():
                if rand() < p:
                    self.eps[e[0],e[1]] = -1*self.eps[e[0],e[1]]
                    self.eps[e[1],e[0]] = self.eps[e[0],e[1]]

    def is_connected(self):
        """
        Returns True if the graph is connected (no disconnected subgraphs), False otherwise.
        """
        return nx.algorithms.is_connected(self)


    def number_of_edges(self):
        """
        Returns the number of edges in the graph.
        """
        return len(self.edges())


    def degree_mean_var(self):
        """
        Returns the mean and variance of the node degrees.
        """
        return mean(self.degree().values()),var(self.degree().values())


    def length_mean_var(self):
        """
        Returns the mean and variance of the edge lengths.
        """
        lsum = 0.0
        l2sum = 0.0
        for e1,e2 in self.edges():
            l = self[e1][e2]['length']
            lsum += l
            l2sum += l*l
        lmean = lsum/self.number_of_edges()
        return lmean,(l2sum/self.number_of_edges() - lmean*lmean)


    def connect_empty(self,N):
        """
        Adds N nodes to the graph, but no edges.  This can be used to clear the graph without deleting
        the object.
        """
        self.remove_nodes_from(self.nodes())
        # re-add desired number of nodes
        self.add_nodes_from(range(0,N))


    def connect_full(self,N,p):
        """
        Each node is connected to every other node; all N nodes have degree N-1.
        """
        self.connect_empty(N)
        self.add_edges_from(nx.random_regular_graph(N-1,N).edges())


    def connect_ring(self,N):
        """
        Each of the N nodes is connected to its 'neighbors' (node N to N-1 and N+1, modulo N).
        The neighbors are only meaningful once a distance embedding is chosen; if the unit circle
        mapping is chosen, this topolgy gives a ring with nearest neighbor connections.
        """
        self.connect_empty(N)
        for n in self.nodes():
            self.add_edge(n,mod(n+1,N))
            self.add_edge(n,mod(n+N-1,N))


    def connect_fixed_degree(self,N,p):
        """
        All nodes have identical degree; they are each randomly connected to p*N other nodes.
        If p > 1 - 1/N, this will return the regular, fully connected graph.'
        """
        self.connect_empty(N)
        d = int(p*N)
        self.add_edges_from(nx.random_regular_graph(d,N).edges())


    def connect_fixed_edges(self,N,p):
        """
        A fixed fraction of the total possible N*(N-1)/2 connections are made. (Not a binomial
        graph!  The number of edges is always p*N(N-1)/2, not just in the N->infinity limit.)
        """
        self.connect_empty(N)
        dN = int(p*N*(N-1)/2)
        self.add_edges_from(nx.gnm_random_graph(N,dN).edges())


    def connect_gnp(self,N,p):
        """
        Erdos-Renyi (Poisson random) graph G(N,p).
        """
        # this is kind of a dumb way to do this
        self.connect_empty(N)
        self.add_edges_from(nx.gnp_random_graph(N,p).edges())


    def connect_gnm(self,N,m):
        """
        Erdos-Renyi (Poisson random) graph G(N,m).
        """
        self.connect_empty(N)
        self.add_edges_from(nx.gnm_random_graph(N,m).edges())


    def connect_barabasi_albert(self,N,m):
        """
        Barabasi-Albert preferential attachment graph with N nodes and m
        edges from each new node to existing nodes.
        """
        if m > N:
            m = N-1
        # again, not the best way to do this
        self.connect_empty(N)
        self.add_edges_from(nx.barabasi_albert_graph(N,m).edges())


    def connect_newman_watts_strogatz(self,N,k,p):
        """
        Newman-Watts-Strogatz graph staring with a k-nearest neighbor ring.  Additional edges are added with
        probability p.
        """
        # ditto
        self.connect_empty(N)
        self.add_edges_from(nx.newman_watts_strogatz_graph(N,k,p).edges())


    def connect_watts_strogatz(self,N,k,p):
        """
        Watts-Strogatz graph starting with a k-nearest neighbor ring.  Edges are then rewired with
        probability p.  This differs from the NWS graph in that (1) it may not be connected and
        (2) the mean degree will be fixed at k.
        """
        self.connect_empty(N)
        self.add_edges_from(nx.watts_strogatz_graph(N,k,p).edges())

    def connect_configuration(self,N,degseq):
        """
        Uses the configuration model to wire up a pulse oscillator network with a given
        input degree sequence. The length of degseq must equal N.  The resulting network
        is pruned of both self loops and parallel (multi) edges.
        """
        assert len(degseq) == N, "ERROR. Each node needs an input degree."
        self.connect_empty(N)
        G = nx.configuration_model(degseq)
        G = nx.Graph(G)
        G.remove_edges_from(G.selfloop_edges())
        self.add_edges_from(G.edges())


    def connect_edgelist(self,N,edgelist):
        """
        Creates the network from an input edgelist.
        """
        self.connect_empty(N)
        self.add_edges_from(edgelist)


    def set_edge_lengths(self,embedding):
        """
        Sets the 'length' attribute of the graph edges according to some physical mapping, stored
        in a DistanceEmbedding object.
        """
        for e1,e2 in self.edges():
            self[e1][e2]['length'] = embedding.distances[e1,e2]
        return


    def euler_integrate(self,p,y0,T,M=10000,fullout=True,stopatsync=False):
        """
        Integrates (using the Euler method) a delayed pulse-oscillator network.
        Currently, only one kind of RHS (Mirollo-Strogatz) is supported.

        INPUT:
            p : array, required
                parameters required for RHS dydt (cannot be none!)

            y0 : array, required
                vector of initial conditions, length equal to number of nodes in network

            T : float, required
                total integration time (integration is performed from 0 to T)

            M : integer, optional
                total number of steps in the integration

            fullout : boolean, optional
                    if fullout == True, an array of size nNodes x M will be returned, giving
                    y(t) at all M simulation steps.  Set to False to return only y(T).

            stopatsync : boolean, optional
                    if stopatsync == True, the integration will terminate as soon as all
                    the nodes are synchronized (i.e., all reset in the same step before
                    pulses are resolved)

        OUTPUT:
            y : amplitudes for each node at all simulation times (output == 'full') or just
                the final integration time (output == 'final')

        OUTPUT:

            from [0,T]
        using M total steps.  Returns the final values for the node amplitudes.
        """
        # lots of massaging to:
        #   1. maintain type consistency with the cython/C call
        #   2. avoid passing too many python objects/functions into the C call
        # booleans to uint8
        fo = uint8(0)
        if fullout is True:
            fo = uint8(1)
        sos = uint8(0)
        if stopatsync is True:
            sos = uint8(1)
        # make sure parameters are floats
        p = p.astype(float64)
        # type and shape of y0
        y0 = asarray(y0).reshape((len(self.nodes()),1)).astype(float64)
        # things to try to avoid lots of python object access
        yth = float64(self.y_th)
        delta = float64(self.delta)
        eps = self.eps.astype(float64)
        #eps = float64(self.eps)
        lengthAdj = zeros((len(self.nodes()),len(self.nodes())),dtype=float64)
        for i in xrange(len(self.nodes())):
            nlist = self.neighbors(i)
            for n in nlist:
                lengthAdj[i,n] = self[i][n]['length']
        # make the call to the integrator
        y,s = eulerint.euler_integrate(lengthAdj,p,y0,yth,delta,eps,T,M,fo,sos)
        return y,s


class DistanceEmbedding(object):
    """
    Calculates node-node distances for a different physical embeddings of a graph.
    Each call to a map function will replace the current distance matrix (distances)
    with a new one, recalculated as described.  If the number of nodes in the graph
    changes, you must re-initialize the embedding as well.

    Methods:

        null_map():
            All distances are equal and of unit magnitude.

        unitcirc_map():
            Equal spacing of the N nodes on a unit circle.

    """
    def __init__(self,N):
        self.N = N
        self.distances = zeros((self.N,self.N))


    def null_map(self):
        """
        All distances are assumed to be equal and unit magnitude (no embedding).
        """
        self.distances = ones((self.N,self.N))
        # zero out self distances
        for i in xrange(self.N):
            self.distances[i,i] = 0.0
        return


    def unitcirc_map(self):
        """
        Assign each node to a position on the unit circle, equally spaced.
        Node 0 is closest to nodes N-1 and 1, 1 closest to 0 and 2, etc.
        """
        dtheta = 2*pi/self.N
        x = cos(arange(self.N)*dtheta)
        y = sin(arange(self.N)*dtheta)
        for i in xrange(self.N):
            for j in xrange(i+1,self.N):
                self.distances[i,j] = sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2)
                self.distances[j,i] = self.distances[i,j]
        return
