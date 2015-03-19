from pydynet.network import DistanceEmbedding
import numpy as np

class TestEmbedding:

    def setup(self):
        self.embedding = DistanceEmbedding(10)

    def test_unitcirc_symmetry(self):
        self.embedding.unitcirc_map()
        assert np.allclose(0,self.embedding.distances - self.embedding.distances.T), 'UNIT CIRCLE: Distance matrix is not symmetric!'


