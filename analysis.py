#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

Created by Kevin Brown on 2015-03-17.
"""
from numpy import log,exp,mean,abs,log2
from numpy import roll,where,histogram,nonzero,delete,zeros_like,array,zeros,newaxis

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
    Computes a dictionary of population codewords from the input spike array.
    Codewords are keyed on word (represented as a string) and entries indicate
    counts.

    INPUT:
        spikes: array-like, required
            should be an N x t array of integers (0 and 1 only)

    OUTPUT:
        codewords: dictionary
            number of times (up to t) that each word appears
    '''
    codewords = {}
    N,t = spikes.shape
    for k in xrange(0,t):
        word = ''.join([str(x) for x in spikes[:,k]])
        if codewords.has_key(word):
            codewords[key] += 1
        else:
            codewords[key] = 1
    return codewords


def convert_to_spikes(y,sorting='lower',thresh=1.0e-06):
    '''
    Accepts an input array y (N x t), representing amplitude or something
    similar, and returns an array of 1's and 0's of the same size giving
    spike locations, defined as times in which the amplitude of a node
    was above or below a user-defined threshold.

    INPUT:
        y: array-like, required
           N x t array; time series for each of the N nodes in the model

        sorting: string, optional
            set to 'upper' to flag times where y(i,t) > thresh, and
            'lower' to flag times where y(i,t) < thresh

        thresh: float, optional
            cutoff for flagging times; spikes can probably be identified as
            times when y(i,t) > thresh or y(i,t) < thresh.  For example, with
            sorting = 'upper', thresh might be 1.0 (depending on the EOM for
            a single node).  If sorting = 'lower', thresh might be 1.0e-06
            (if we are looking for times when the node has reset to a value
            of zero).

    OUTPUT:
        s : array-like
            N x t array of 1's and 0's, with 'interesting' times in y having
            the value 1
    '''
    s = zeros_like(y,dtype=int)
    if sorting is 'upper':
        s[y >= thresh] = 1
    else:
        s[y <= thresh] = 1
    return s


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


def random_lz_complexity(n,p=0.5):
    '''
    Computes the expected Lev-Zimpel complexity for a random sequence of length
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
    Lev-Zimpel complexity as described in Kaspar and Schuster, Phys. Rev. A.
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
            # determine probability of generating a 1
            p = (sum(spike_array[i,:]) + 1.0)/(T + 2.0)
            # convert the list of spikes to a string
            s = ''.join([str(x) for x in spike_array[i,:]])
            # compute normalized LZ complexity
            c[i] = 1.0*lz_complexity(s)/random_lz_complexity(T,p)
    return c
