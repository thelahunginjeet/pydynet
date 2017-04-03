#!/usr/bin/env python
# encoding: utf-8
from numpy.random import randint

def randchoice(l):
    """
    Returns an element from list l chosen at random.
    """
    return l[randint(len(l))]

def randbit(size=None):
    '''
    Generates an array of shape size of random {0,1} bits.
    '''
    return randint(2,size=size)

def randspin(size=None):
    '''
    Generates an array of shape size of random {-1,1} spin variables.
    '''
    return 2*randbit(size=size) - 1
