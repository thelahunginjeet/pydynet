#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

Created by Kevin Brown on 2015-03-17.
"""
from numpy import log,histogram,nonzero,delete,zeros_like,array,exp,mean,abs

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


def convert_to_spikes(y,sorting='upper',thresh=1.0):
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
            the value 1.0
    '''
    s = zeros_like(y)
    if sorting is 'upper':
        s[y >= thresh] = 1
    else:
        s[y <= thresh] = 1
    return s


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
