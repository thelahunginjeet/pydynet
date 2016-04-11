"""
plotting.py

This module has wrappers for making figures that we find we often need.

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

"""

import pylab
import networkx as nx
from numpy import ceil,exp

def plot_spike_raster(spike_raster):
    '''
    Accepts an input integer array of 1's and 0's (a 1 denoting a spike fired
    in that time bin) and uses imshow, with reasonable options to produce a
    b/w image.  Spikes are shown in black.
    '''
    fig = pylab.figure()
    # the 1-spike_raster converts zeros (no spikes) to white
    pylab.imshow(1-spike_raster,interpolation='none',aspect='auto',cmap='gray',figure=fig)
    fig.axes[0].set_xticks([])
    fig.axes[0].set_yticks([])
    return fig

def plot_network_ring(G,defcolor='k',ncData=None,layout='radial'):
    '''
    Draws an input graph G by arranging the nodes on a ring. Returns a pylab
    figure object for further manipulation or saving/display.

    INPUT
    ------
    G     : networkx graph, required

    defcolor : string, optional
        default color for nodes

    ncData  : dict, optional
        custom color information to decorate nodes; should be a dictionary
        keyed on nodes, with valid color information (strings, floats etc.).
        Only nodes which should have the nondefault color need to be present.

    layout  : string, optional
        suggested layouts are 'circo','twopi',and 'radial'
    '''
    # construct node colors
    nList = G.nodes()
    nc = [defcolor for x in nList]
    if ncData is not None:
        for i in xrannge(0,len(nList)):
            if ncData.hask_key(nList[i]):
                nc[i] = ncData[nList[i]]
    # now make the plot
    fig = pylab.figure(figsize=(8,8))
    if layout == 'radial':
        pos = nx.circular_layout(G)
    else:
        pos = nx.nx_pydot.graphviz_layout(G,prog=layout)
    # heuristic attempt to automatically determine node size
    n = len(G.nodes())
    ns = max(ceil(300*exp(-0.00856*(n - 10))),25)
    ns = int(ns/2.0)
    # make the figure
    nx.draw_networkx(G,pos=pos,node_color=nc,node_size=ns,with_labels=False,figure=fig)
    fig.axes[0].set_axis_off()
    return fig
