#!/usr/bin/env python
# encoding: utf-8
"""
analysis.py

@author: Kevin S. Brown (University of Connecticut), Ann M. Hermundstad (UPenn)

Created by Kevin Brown on 2015-03-17.
"""
from numpy import log,histogram,nonzero,delete

def entropy(x,bins=10,est='ML'):
    """
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
    """
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