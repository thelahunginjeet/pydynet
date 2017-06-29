#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

Created by Kevin Brown on 2015-03-17.
"""
from numpy import log,exp,mean,abs,log2,sqrt,dot,power,log10,logspace,median
from numpy import roll,where,histogram,nonzero,delete,zeros_like,array,zeros,newaxis,array_split
from numpy import append,insert
import networkx as nx

def phi_of_t(y,group=None):
    '''
    For y(i,t), computes:

        phi(t) = |<exp(i*y(t))>|

    where the average is computed only over the rows indicated in group.  If
    group is None, the average is computed over the entire array.

    INPUT:
        y: array-like, required
            y should be an N x t array of floats/doubles.

        group : list, optional
            group should be a list of rows of y to include in the calculation.
            if group is None, all rows are used.

    OUTPUT:
        phi : array
            phi will be an array of size 1 x t, representing the synchonization
            index, as a function of time, computed over the input row group
    '''
    if group is None:
        group = range(0,y.shape[0])
    # set up the array exp(i*y)
    eiy = exp(1j*y)
    # do the averaging
    phi = eiy[group,:].mean(axis=0)
    # now the modulus
    return abs(phi)


def codeword_dictionary(spikes):
    '''
    Computes a dictionary of population codewords from an input spike array.
    Codewords are keyed on word (represented as a string) and entries indicate
    counts.

    INPUT:
        spikes: array-like, required
            should be an N x t array of integers (0 and 1 only)

    OUTPUT:
        codewords: dictionary
            number of times (up to t) that each word appears

        codeseq: list
            order of appearance of codewords; each codeword is assigned an
            arbitrary number between 0 and N(codewords)-1

        codetonum: dictionary
            codeword to numerical mapping in codeseq
    '''
    codeseq = []
    codetonum = {}
    codewords = {}
    N,t = spikes.shape
    current = 0
    for k in xrange(0,t):
        word = ''.join([str(x) for x in spikes[:,k]])
        if codewords.has_key(word):
            codewords[word] += 1
        else:
            codewords[word] = 1
            codetonum[word] = current
            current += 1
        codeseq.append(codetonum[word])
    return codewords,codetonum,codeseq


def codeword_raster(codewords):
    '''
    Uses a codeword dictionary to produce a codeword raster (codewords x spikes)
    and a numeric array of codeword frequencies.
    '''
    coderast = []
    codenum = []
    for k in codewords:
        coderast.append([int(c) for c in k])
        codenum.append(1.0*codewords[k]/sum(codewords.values()))
    coderast = array(coderast)
    codenum = array(codenum)
    return codenum,coderast

def bin_spikes(spike_array,b=10):
    '''
    Accepts an input integer array of zeros and ones and bins samples along
    the column index.  Binned samples are computed according to
        value(bin) = max(bin).
    '''
    n,t = spike_array.shape
    binned_array = zeros((n,t/b),dtype=int)
    binstart = 0
    while t >= binstart + b:
        binvals = spike_array[:,binstart:binstart+b].max(axis=1)
        binned_array[:,binstart:binstart+b] = binvals[:,newaxis]
        binstart += b
    return binned_array


def isi_stats(spike_array):
    '''
    Accepts an N x t input array of 1's and 0's, with a 1 indicating a spike
    occurred in that time bin and returns the mean and variance of the interspike
    intervals.

    INPUT:
        spike_array : array, required
            spike_array should only contain 1's and 0's; amplitude/phase arrays
            should be pre-converted via convert_to_spikes

    OUTPUT:
        isi_mean : array (N elements)
            ISI mean

        isi_var : array (N elements)
            ISI variance
    '''
    N = spike_array.shape[0]
    isi_mean = zeros(N)
    isi_var = zeros(N)
    for k in xrange(0,N):
        spike_loc = where(spike_array[k,:] == 1)[0]
        isi_array = spike_loc - roll(spike_loc,1)
        isi_mean[k] = isi_array[1:].mean()
        isi_var[k] = isi_array[1:].var()
    return isi_mean,isi_var


def discrete_entropy(x,est='ML'):
    '''
    Computes the entropy of discrete (integer) data.

    INPUT:
        x: array-like, required
            data for which to compute H[x]

        est: string, optional
            estimator. current options arange

                'ML' : maximum-likelihood (plugin)
                'MM' : Miller-Maddow corrected

    OUTPUT:
        H[x]: entropy of x, measured in nats
    '''
    # do the frequency counting
    counts = {}
    for xi in x:
        if counts.has_key(xi):
            counts[xi] += 1
        else:
            counts[xi] = 1
    sumpofx = 1.0*sum(counts.values())
    pofx = array(counts.values())/sumpofx
    H_ML = -1*(pofx*log(pofx)).sum()

    if est == 'ML':
        H = H_ML

    if est == 'MM':
        # nonzero bins have already been removed from pofx
        H = H_ML + (len(pofx) - 1.0)/(2.0*len(x))
    return H


def entropy(x,bins=10,est='ML'):
    '''
    Computes the entropy of a continuous (unbinned) set of data x.

    INPUT:
        x: array-like, required
           data for which to compute H[x]

        bins: integer, optional
           number of bins

        est: string, optional
           estimator.  current options are:

               'ML':  maximum-likelihood (plugin)
               'MM':  Miller-Maddow corrected
               'JK':  Jackknifed estimate (can be slow!)

    OUTPUT:
        H[x] : entropy of x, measured in nats
    '''
    cx = histogram(x,bins)[0]
    pofx = (1.0*cx)/cx.sum()
    # remove zero bins to avoid numerical problems
    pofx = pofx[nonzero(pofx)]
    H_ML = -1*(pofx*log(pofx)).sum()

    if est == 'ML':
        H = H_ML

    if est == 'MM':
        # nonzero bins have already been removed from pofx
        H = H_ML + (len(pofx) - 1.0)/(2.0*len(x))

    if est == 'JK':
        Sc = 0
        for i in xrange(0,len(x)):
            newx = delete(x,i)
            Sc += entropy(newx,bins,'ML')
        H_JK = len(x)*H_ML - ((len(x) - 1.0)/len(x))*Sc
        H = H_JK

    return H


def codeword_complexity(spike_array,norm=True):
    '''
    Computes the Lempel-Ziv complexity for a series of codewords.  If norm is
    True, the normalized lz_complexity is returned.  Also returns the number
    of unique codewords.
    '''
    N,t = spike_array.shape
    # find and count the codewords
    codewords,codetonum,codeseq = codeword_dictionary(spike_array)
    nunique = len(codetonum.keys())
    # compute the non-normalized LZ complexity
    lzc = lz_complexity(codeseq)
    # normalize if desired
    if norm is True:
        f = 1.0*array(codewords.values())/t
        # source entropy
        h = -sum(f*log2(f))
        # length term
        bn = t/log2(t)
        # normalize
        lzc = lzc/(h*bn)
    return lzc,nunique


def random_lz_complexity(n,p=0.5):
    '''
    Computes the expected Lempel-Ziv complexity for a random sequence of length
    n and expected probability of generating a 1 = p.  Useful for normalizing
    the raw lz_complexity.  This function will behave poorly if p is identically
    0 or 1.  Therefore, it would be best to estimate p from real (finite length)
    strings using pseudocounts.

    INPUT:
        n : int, required
          length of the random sequence

        p : float, optional
          probability of seeing a 1 in the sequence
    '''
    # source entropy
    h = -p*log2(p) - (1-p)*log2(1-p)
    # expected LZ complexity of binary representations of real numbers
    bn = n/log2(n)
    return h*bn



def lz_complexity(s):
    '''
    Lempel-Ziv complexity as described in Kaspar and Schuster, Phys. Rev. A.
    The input iterable (see below) does not have to be binary (2-element), but
    most applications of LZ complexity have used strings of 0s and 1s.

    INPUT:
        s : string, list, or tuple, required
          sequence to calculate complexity for

    '''
    i, k, l = 0, 1, 1
    k_max = 1
    n = len(s)-1
    lzc = 1
    while True:
        if s[i+k-1] == s[l+k-1]:
            k += 1
            if l + k >= n - 1:
                lzc += 1
                break
        else:
            if k > k_max:
               k_max = k
            i += 1
            if i == l:
                lzc += 1
                l += k_max
                if l + 1 > n:
                    break
                else:
                    i = 0
                    k = 1
                    k_max = 1
            else:
                k = 1
    return lzc


def complexity(spike_array,method='lz_norm'):
    '''
    Complexity measure for each node's spiking pattern.  Could dispatch to a
    variety of measures.  Returns an array of length equal to spike_array.shape[0].
    '''
    N,T = spike_array.shape
    c = zeros(N)
    if method == 'lz_norm':
        for i in xrange(0,N):
            # spike string
            s = ''.join([str(x) for x in spike_array[i,:]])
            # probability of generating a 1
            p = (sum(spike_array[i,:]) + 1.0)/(T + 2.0)
            # compute normalized LZ complexity
            c[i] = 1.0*lz_complexity(s)/random_lz_complexity(T,p)
    if method == 'lz':
        for i in xrange(0,N):
            # spike string
            s = ''.join([str(x) for x in spike_array[i,:]])
            # non-normalized lz complexity
            c[i] = lz_complexity(s)
    return c


def rolling_complexity(s,block_length=2000,block_shift=1000):
    '''
    Computes rolling complexity (and IPR) for a spike raster s.  Calculations
    are performed using overlapping blocks of size block_length that shift
    in steps of block_shift.

    INPUT:
        s: integer array, required
            binary array of nodes x times

        block_length : integer, optional
            block size for LZC calculation

        block_shift : integer, optional
            number of steps to advance the blocks

    OUTPUT:
        t : array
            array of times at the end of blocks

        nlzc : array
            nLZC for each block

        ipr : array
            Inverse Participation Ratio for each block
    '''
    nlzc = []
    ipr = []
    t = []
    t_start = block_length
    next_start = 0
    while True:
        x = s[:,next_start:next_start + block_length]
        lzc = complexity(x)
        nlzc.append(median(lzc))
        ipr.append(inv_part_ratio(lzc))
        t.append(next_start+block_length)
        next_start += block_shift
        if next_start + block_length > s.shape[1]:
            break
    return array(t),array(nlzc),array(ipr)


def node_assortativity(net,attribute,jackknife=True,atype='numeric'):
    '''
    Computes the assortativity coefficient and optional sampling error (via the
    jackknife) for the desired attribute over the network net.  In addition, this
    only works as expected for unweighted, undirected graphs.

    This function assumes the input nodes are not already decorated with the attribute;
    this will almost always be the case when the attribute arises as a post-simulation
    calculation on the dynamics of the network.

    INPUT:
        net: PulseOscillatorNetwork (or networkx graph), required
            input network

        attribute : dictionary, required
            key/value pairs for the attribute; keys should be valid node names
            in the network net

        jackknife : bool, optional
            set to True to compute the expected sampling variance

        atype : string, optional
            set to 'numeric' for integer point attributes and
            'categorial' for categorical attributes

    OUTPUT:
        r : float
            numerical attribute assortativity coefficient (-1 < r <= 1)

        sigmar : float, optional
            jackknife standard deviation of r
    '''
    # set the correct assortativity function
    if atype is 'numeric':
        afunc = nx.numeric_assortativity_coefficient
    else:
        afunc = nx.attribute_assortativity_coefficient
    # create a new graph
    G = nx.Graph()
    # add nodes from the network, with attributes
    for n in net.nodes():
        G.add_node(n,value=attribute[n])
    G.add_edges_from(net.edges())
    r = afunc(G,'value')
    if jackknife:
        sigmarsq = 0.0
        # remove one edge at a time, recompute, then add it back
        for e in G.edges():
            G.remove_edge(e[0],e[1])
            sigmarsq += (afunc(G,'value') - r)**2
            G.add_edge(e[0],e[1])
        return r,sqrt(sigmarsq/len(G.edges()))
    else:
        return r


def sparsity(v):
    '''
    Computes sparsity (Vinje and Gallant, Science 287(18) 2000) for a vector v.
    Defined as:
            S = (1 - <v>^2/<v^2>)(1 - 1/n)^{-1}
    where n = len(v).
    '''
    n = len(v)
    return (1.0 - (v.mean()**2)/((v**2).mean()))/(1.0 - 1.0/n)


def inv_part_ratio(v):
    '''
    Computes the inverse participation ratio for vector v; v is normalized
    before calculation.
    '''
    vnorm = sqrt(dot(v,v))
    v2 = power(v/vnorm,2)
    v4 = power(v/vnorm,4)
    return power(v2.sum(),2)/v4.sum()

def fano_factor_tc(s):
    '''
    Computes the Fano factor for an input spike train s as a function of the
    block length.  The block lengths and Fano factor for those block lengths
    are both returned.  The Fano factor for block length T is defined as:

                    F(T) = var(N_i(T))/<N_i(T)>
    '''
    # this is the number of pieces to cut
    nchunks = (len(s)/logspace(0.0,log10(1.0*len(s)/10),num=20)).astype(int)
    fano = zeros(len(nchunks))
    for i in range(len(nchunks)):
        Ni = array([sum(x) for x in array_split(s,nchunks[i])])
        fano[i] = Ni.var()/Ni.mean()
    return 1.0*len(s)/nchunks,fano


def allan_factor_tc(s):
    '''
    Computes the Allan factor for an input spike train s as a function of the
    block length.  The block lengths and Allan factor for those block lengths
    are both returned.  The Allan factor for block length T is defined as:

                    A(T) = <(N_(i+1)(T) - N_i(T))^2>/2<N_i(T)>
                         = 2F(T) - F(2T)
    '''
    nchunks = (len(s)/logspace(0.0,log10(1.0*len(s)/10),num=20)).astype(int)
    allan = zeros(len(nchunks))
    for i in range(len(nchunks)):
        Ni = array([sum(x) for x in array_split(s,nchunks[i])])
        anum = array([(Ni[j+1]-Ni[j])**2 for j in range(0,len(Ni)-1)])
        allan[i] = anum.mean()/2*Ni.mean()
    return 1.0*len(s)/nchunks,allan


def k_teager(x,k):
    '''
    Computes the k-Teager energy operator for input vector x.  The k-Teager
    energy is defined as:

            E[n] = max(x[n]**2 - x[n-k]*x[n+k],0)
    '''
    # pad at the beginning and end of the array
    x = append(x,[x[-1]]*k)
    x = insert(x,0,[x[0]]*k)
    # compute
    E = x**2 - roll(x,k)*roll(x,-k)
    E[E < 0] = 0.0
    return E
