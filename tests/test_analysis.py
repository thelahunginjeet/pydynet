from pydynet.analysis import *
import numpy as np

class TestAnalysis:

    def setup(self):
        pass

    def test_entropy_ML(self):
        x = np.random.randn(1000)
        H = entropy(x,10,'ML')
        assert H > 0, 'ML ENTROPY: Entropy is negative!'

    def test_entropy_MM(self):
        x = np.random.randn(1000)
        H = entropy(x,10,'MM')
        assert H > 0, 'MM ENTROPY: Entropy is negative!'

    def test_entropy_JK(self):
        x = np.random.randn(1000)
        H = entropy(x,10,'JK')
        assert H > 0, 'JK ENTROPY: Entropy is negative!'