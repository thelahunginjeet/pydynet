"""
plotting.py

This module has wrappers for making figures that we find we often need.

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

"""

import pylab
import networkx as nx
from numpy import ceil,exp,array,where,arange,hstack,ones


def plot_spike_raster(spike_array,figtype='image',msize='9'):
    '''
    Accepts an input integer array of 1's and 0's (a 1 denoting a spike fired
    in that time bin) and uses either:
        (1) figtype:'image' : imshow, with reasonable options to produce a b/w image.
            Spikes are shown in black.
        (2) figtype:'dots' : plot, in which each spike is represented as a '.'  size
            of the dot is controlled via msize.
    '''
    fig = pylab.figure()
    if figtype is 'image':
        # the 1-spike_raster converts zeros (no spikes) to white
        pylab.imshow(1-spike_array,interpolation='none',aspect='auto',cmap='gray',figure=fig)
    else:
        # spikes are present where there is a '1' in the spike_array
        N = spike_array.shape[0]
        t = arange(spike_array.shape[1])
        for i in range(N):
            spikeloc = where(spike_array[i,:] == 1)[0]
            pylab.plot(t[spikeloc],i*spike_array[i,spikeloc],'k.',markersize=msize)
            pylab.xlim([-1,len(t)+1])
            pylab.ylim([-0.5,N-0.5])
    fig.axes[0].set_xticks([])
    fig.axes[0].set_yticks([])
    return fig


def construct_ncdata(G,nlists,rgbcolors,defcolor=[0.,0.,0.]):
    '''
    Constructs a dictionary of node colors for use in plot_network_ring.
    listsofnodes.  rgbcolors should be a
    list of 3-lists, with row dimension equal to the length of args.
    rgbcolors[i,:] will be assigned to args[i].  Any nodes not
    specified in the list of args will receive the default color.

    Example:

        construct_ncdata(G,[[1,2,3],[4,5,6]],[[1.,0.,0.],[0.,1.,0.]])

    Will set up node coloring so that [1,2,3] are red, [4,5,6] are blue, and
    all other nodes are black.
    '''
    ncData = {}.fromkeys(G.nodes())
    for k in ncData:
        ncData[k] = defcolor
    for i in range(len(nlists)):
        for n in nlists[i]:
            ncData[n] = rgbcolors[i]
    return ncData


def construct_ecdata(G,nodelist,lightcolor=0.25):
    '''
    Constructs a list of edgecolors for use in plot_network_ring.  All edges
    connected to the nodes in nodelist will be given a color of 1.0 (what color
    that is depends on the edge colormap in plot_network_ring), and the other
    edges will have a color of 0.5.  This is most useful for "graying out"
    edges not connected to the nodes in nodelist.
    '''
    edgelist = G.edges()
    ecData = lightcolor*ones(len(edgelist))
    for i in range(len(ecData)):
        e = edgelist[i]
        for n in nodelist:
            if e[0] == n or e[1] == n:
                ecData[i] = 1.0
    return ecData


def plot_network_ring(G,nodeline=False,defcolor='k',ncData=None,ecData=None,layout='radial',cmap=pylab.cm.plasma,edgecmap=pylab.cm.Greys):
    '''
    Draws an input graph G by arranging the nodes on a ring. Returns a pylab
    figure object for further manipulation or saving/display.

    INPUT
    ------
    G     : networkx graph, required

    nodeline : boolean, optional
        draw node outlines or suppress them?

    defcolor : string, optional
        default color for nodes

    ncData  : dict, optional
        custom color information to decorate nodes; should be a dictionary
        keyed on nodes, with valid color information (strings, floats etc.).
        Nodes not in the dictionary receive the default color.

    ecData : list, optional
        custom color information to decorate edges; should be a list/array
        of color information ordered the same way as G.edges().  If set to
        none, all edges are drawn in black.

    layout  : string, optional
        suggested layouts are 'circo','twopi',and 'radial'
    '''
    # construct node colors
    nList = G.nodes()
    nc = [defcolor for x in nList]
    if ncData is not None:
        for i in xrange(0,len(nList)):
            if ncData.has_key(nList[i]):
                nc[i] = ncData[nList[i]]
    # now make the plot
    fig = pylab.figure()
    if layout == 'radial':
        pos = nx.circular_layout(G)
    else:
        pos = nx.nx_pydot.graphviz_layout(G,prog=layout)
    # heuristic attempt to automatically determine node size
    n = len(G.nodes())
    ns = max(ceil(300*exp(-0.00856*(n - 10))),25)
    ns = int(ns/2.0)
    # edge colors
    if ecData is None:
        ecData = 'k'
    # make the figure
    if nodeline is True:
        nx.draw_networkx(G,pos=pos,node_color=nc,node_size=ns,with_labels=False,figure=fig,cmap=cmap,edge_color=ecData,edge_cmap=edgecmap,edge_vmin=0.0,edge_vmax=1.0)
    else:
        nodes = nx.draw_networkx_nodes(G,pos=pos,node_color=nc,node_size=ns,with_labels=False,figure=fig,cmap=cmap)
        nodes.set_edgecolor('none')
        nx.draw_networkx_edges(G,pos=pos,edge_color=ecData,edge_cmap=edgecmap,edge_vmin=0.0,edge_vmax=1.0,figure=fig)
    ax = pylab.gca()
    ax.set_aspect('equal')
    fig.axes[0].set_axis_off()
    return fig
