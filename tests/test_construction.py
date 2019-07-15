from __future__ import division
from builtins import object
from past.utils import old_div
from pydynet.network import PulseOscillatorNetwork
import numpy as np

'''Test for other values of N.'''


class TestConstruction(object):

    def setup(self):
        self.N = 10

    def test_connect_empty(self):
        net = PulseOscillatorNetwork(self.N,topology='empty')
        assert len(net.nodes()) == self.N, 'EMPTY: Wrong number of nodes!'
        assert len(net.edges()) == 0, 'EMPTY: Edges present!'

    def test_connect_full(self):
        net = PulseOscillatorNetwork(self.N,topology='full')
        assert len(net.nodes()) == self.N, 'FULLY CONNECTED: Wrong number of nodes!'
        assert len(net.edges()) == old_div(self.N*(self.N-1),2), 'FULLY CONNECTED: Wrong number of edges!'
        assert np.allclose(self.N-1,list(net.degree().values())), 'FULLY CONNECTED: Some nodes have the wrong degree!'

    def test_connect_ring(self):
        net = PulseOscillatorNetwork(self.N,topology='ring')
        assert len(net.nodes()) == self.N, 'RING: Wrong number of nodes!'
        assert len(net.edges()) == self.N, 'RING: Wrong number of edges!'
        assert np.allclose(2,list(net.degree().values())), 'RING: Some nodes have the wrong degree!'


    def test_connect_fixed_degree(self):
        p = np.random.uniform(high=1-1./self.N)
        net = PulseOscillatorNetwork(self.N,p=p,topology='fixed degree')
        assert len(net.nodes()) == self.N, 'FIXED DEGREE: Wrong number of nodes!'
        assert len(net.edges()) == old_div((int(p*self.N)*self.N),2), 'FIXED DEGREE: Wrong number of edges!'
        assert np.allclose(int(p*self.N),list(net.degree().values())), 'FIXED DEGREE: Some nodes have the wrong degree!'

    def test_connect_fixed_edges(self):
        p = np.random.rand()
        net = PulseOscillatorNetwork(self.N,p=p,topology='fixed edges')
        assert len(net.nodes()) == self.N, 'FIXED EDGES: Wrong number of nodes!'
        assert len(net.edges()) == int(old_div(p*self.N*(self.N-1),2)), 'FIXED EDGES: Wrong number of edges!'

