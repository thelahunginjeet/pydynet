from pydynet import rewiring
from pydynet.network import PulseOscillatorNetwork

class TestRewiring:

    def setup(self):
        pass

    def test_add_edge(self):
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        nEdges = net.number_of_edges()
        edges = net.edges()
        flag = rewiring.add_random_edge(net)
        if flag:
            assert net.number_of_edges() == nEdges + 1, "ADD: Edge should have been added!"
        else:
            assert net.edges() == edges, "ADD: Addition should have failed!"

    def test_remove_edge(self):
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        nEdges = net.number_of_edges() 
        edges = net.edges()
        flag = rewiring.remove_random_edge(net)
        if flag:
            assert net.number_of_edges() == nEdges - 1, "REMOVE: Edge should have been removed!"
        else:
            assert net.edges() == edges, "REMOVE: Edges should not have changed!"

    def test_move_edge(self):
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        nEdges = net.number_of_edges() 
        edges = net.edges()
        flag = rewiring.move_random_edge(net)
        if flag:
            assert net.number_of_edges()== nEdges, "MOVE: Number of edges should not have changed!"
        else:
            assert net.edges() == edges, "MOVE: Edges should not have changed!"

    def test_swap_edges(self):
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        nEdges = net.number_of_edges() 
        edges = net.edges()
        flag = rewiring.swap_random_edges(net)
        if flag:
            assert net.number_of_edges() == nEdges, "SWAP: Number of edges should not have changed!"
        else:
            assert net.edges() == edges, "SWAP: Edges should not have changed!"

    def test_move_degree_dist(self):
        from numpy import sort
        accept = 0.0
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        nEdges = net.number_of_edges()
        for i in xrange(0,100):
            accept += rewiring.move_random_edge(net)
        assert net.number_of_edges() == nEdges, "MOVE DEGREE: Number of edges should not have changed!"
        
    def test_swap_node_degrees(self):
        accept = 0.0
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        deg = net.degree()
        for i in xrange(0,100):
            accept += rewiring.swap_random_edges(net)
        assert accept > 0.0, "SWAP DEGREE: No swaps were performed!"
        for k in net.degree():
            assert net.degree()[k] == deg[k], "SWAP DEGREE: Node degrees should not have changed!"

    def test_graph_perturb(self):
        net = PulseOscillatorNetwork(10,0.3,'fixed degree')
        pertNet,R = rewiring.perturb_graph(net)
        assert R > 0.0, "PERTURB: No graph operations performed!"
        assert net != pertNet, "PERTURB: Input graph copy failed!"