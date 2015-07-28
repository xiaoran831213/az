import os
import sys
import time
import os.path as pt
import numpy as np
import theano
import theano.tensor as T
from theano import function as F
from theano.tensor.shared_randomstreams import RandomStreams
from theano import pp
import hlp
from hlp import S
import trainer
from trainer import Trainer
import lyr
from lyr import Lyr
import pdb

def get_data():
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    x = x.reshape(x.shape[0], -1)
    d = (x.shape[1], x.shape[1]/2)
    x = (x - x.min()) / (x.max() - x.min())
    return x

def get_das(dims, data = None):
    rnd = np.random.RandomState(150)
    if data is not None:
        dims = [data.shape[1]] + dims
    ZDT = zip(dims[:-1], dims[1:], xrange(len(dims)))
    ecs = [Lyr((i, j), np_rnd = rnd, tag = 'E{}'.format(t)) for i, j, t in ZDT]
    dcs = [Lyr((j, i), np_rnd = rnd, tag = 'D{}'.format(t)) for i, j, t in ZDT]

    das = zip(ecs, dcs)
    ## constraint weight terms
    for ec, dc in das:
        dc.w(ec.w, T.transpose)
    return das

def pre_train(data, das, nep = 30):
    x = data
    for ec, dc in das:
        dc.x(ec.y)
        tr = Trainer(
            ec.x, dc.y,
            src = x, xpt = x,
            lrt = 0.01, mmt = 0.02,
            call_wreg = trainer.wreg_l2(.1))
        tr.tune(nep)
        ec.x(x)
        x = ec.y().eval()
    del x

def fine_tune(data, das, nep = 30):
    x = data

    ## re-wire encoders and decoders
    ecs, dcs = zip(*das)
    sda = list(ecs) + list(reversed(dcs))
    for i, j in zip(sda[:-1], sda[1:]):
        j.x(i.y) # lower output -> higher input

    tr = Trainer(sda[0].x, sda[-1].y, src = data, xpt = data, lrt = 0.001)
    tr.tune(nep)
    return tr
    
