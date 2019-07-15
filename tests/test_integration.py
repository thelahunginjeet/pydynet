from builtins import object
from pydynet.network import PulseOscillatorNetwork
from pydynet.network import DistanceEmbedding
from pydynet.network import dydtMS
import numpy as np

class TestIntegration(object):

    def setup(self):
        self.net = PulseOscillatorNetwork(10,topology='ring')
        embed = DistanceEmbedding(len(self.net.nodes()))
        embed.unitcirc_map()
        self.net.set_edge_lengths(embed)


    def test_integrate_nodelay(self):
        self.net.eps = 0.1
        self.net.delta = 0.0
        p = np.array([1,2])
        y0 = np.zeros((len(self.net.nodes()),1))
        y = self.net.euler_integrate(p,y0,100,M=10000,fullout=False)


    def test_integrate_delay(self):
        self.net.eps = 0.1
        self.net.delta = 1.0
        p = np.array([1,2])
        y0 = np.zeros((len(self.net.nodes()),1))
        y = self.net.euler_integrate(p,y0,100,M=100000,fullout=False)


    def test_integrate_fullout(self):
        self.net.eps = 0.1
        self.net.delta = 0.0
        p = np.array([1,2])
        y0 = np.zeros((len(self.net.nodes()),1))
        y = self.net.euler_integrate(p,y0,10,M=10000,fullout=True)
        assert y.shape == (len(self.net.nodes()),10000 + 1)


    def test_integrate_truncated(self):
        self.net.eps = 0.1
        self.net.delta = 0.0
        p = np.array([1,2])
        y0 = np.zeros((len(self.net.nodes()),1))
        y = self.net.euler_integrate(p,y0,10,M=10000,fullout=False)
        assert y.shape == (len(self.net.nodes()),1)
