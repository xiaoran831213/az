import os
import sys
import time
import os.path as pt
import numpy as np
import lyr
from lyr import Lyr
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
import hlp
import pdb

def test_da1(ec = None, dc = None, x = None):

    if x is None:
        x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
        x = x.reshape(x.shape[0], -1)
        d = (x.shape[1], x.shape[1]/2)
        x = (x - x.min()) / (x.max() - x.min())

    rnd = np.random.RandomState(120)
    if ec is None:
        ec = Lyr(d = (d[0]/1, d[0]/2), np_rnd = rnd, x = x)
    if dc is None:
        dc = Lyr(d = (d[0]/2, d[0]/1), np_rnd = rnd, x = ec.y)
        dc.w(ec.w, T.transpose)

    return ec, dc

if __name__ == '__main__':
    theano.config.floatX = 'float32'
    pass
