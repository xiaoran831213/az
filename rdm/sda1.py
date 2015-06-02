import os
import sys
import time

import numpy as np

import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams

#from logistic_sgd import LogisticRegression, load_data
import da1
FT = theano.config.floatX

# start-snippet-1
class SDA(list):
    """Stacked denoising auto-encoder class (SdA)

    """

    def __init__(self, n_vis, np_rng, th_rng=None):
        """ This class is made to support a variable number of layers. """

        self.l_da = []
        self.parm = []

        if not th_rng:
            th_rng = RandomStreams(np_rng.randint(2 ** 30))

        self.n_vis = n_vis                # number of visible feature
        self.np_rng = np_rng
        self.th_rng = th_rng
        

    def add_DA(self, n_hid):
        if len(self.l_da):
            n_vis = self.l_da[-1].n_vis
        else:
            n_vis = self.n_vis

        da = da1.DA(
            np_rng = self.np_rng,
            th_rng = self.th_rng,
            n_vis = n_vis,
            n_hid = n_hid)
        
        self.l_da.append(da)
        self.parm.append(da.parm)
        
    def __get_stack__ (self, lyr = None, dp_enc = None, dp_dec = None):
        if not lyr:
            lyr = 0
            dp_enc = -1
            dp_dec = -1
        if lyr < 0:
            lyr = len(self.l_da) - lyr
        if dp_enc == None:
            dp_enc = 1
        if dp_enc < 0:
            dp_enc = sys.maxint

        if dp_dec == None:
            dp_dec = 1
        if dp_dec < 0:
            dp_dec = sys.maxint

        l_h = min(lyr + dp_enc, len(self.l_da))  # last hidden layer + 1
        l_z = max(l_h - dp_dec, 0)               # last re-con layer
        enc = self.l_da[lyr:l_h]
        dec = self.l_da[l_z:l_h]
        dec.reverse()
        return enc, dec
        
    def __get_parms__(self, lyr = None, dp_enc = None, dp_dec = None):
        enc, dec = self.__get_stack__(lyr, dp_enc, dp_dec)
        parms = []
        for da in enc:
            parms.append(da.t_w)
            parms.append(da.t_b)
        for da in dec:
            parms.append(da.t_w_prime)
            parms.append(da.t_b_prime)
        return parms

    def t_pipe(self, t_x, lyr = None, dp_enc = None, dp_dec = None):
        enc, dec = self.__get_stack__(lyr, dp_enc, dp_dec)

        for da in enc:
            t_x = da.t_encode(t_x)
        for da in dec:
            t_x = da.t_decode(t_x)
        return t_x

    def t_encode(self, t_x, lyr = None, dep = None):
        return self.t_pipe(t_x, lyr, dp_enc = dep, dp_dec = 0)

    def t_decode(self, t_h, lyr = None, dep = None):
        return self.t_pipe(t_x, lyr, dp_enc = 0, dp_dec = dep)

    def t_corrupt(self, t_x, lvl):
        return self.th_rng.binomial(
            size = t_x.shape, n = 1, p = 1 - lvl,
            dtype = T.config.floatX) * t_x

    def f_trainer(self, t_data, lv_corrupt = 0.2, t_rate = 0.1,
                  lyr = None, dp_enc = None, dp_dec = None):
        x = T.matrix('x')
        q = self.t_corrupt(t_x, lv_corrupt)
        z = self.t_pipe(x, lyr, dp_enc, dp_dec)

        L = - T.sum(x * T.log(z) + (1 - x) * T.log(1 - z), axis=1)
        cost = T.mean(L)

        dist = T.mean(T.sqrt(T.sum((x - z) ** 2, axis = 1)))

        parm = self.__get_parms__(lyr, dp_enc, dp_dec)

        grad = T.grad(cost, parm)

        diff = [(p, p - t_rate * g) for p, g in zip(parm, grad)]

        t_fr = T.iscalar()
        t_to = T.iscalar()
        return theano.function(
            [t_fr, t_to],
            [cost, dist],
            updates = diff,
            givens = {x : t_data[t_fr:t_to]},
            name = "DA_trainer")

def test_sda():
    np_rng = np.random.RandomState(123)
    th_rng = RandomStreams(np_rng.randint(2 ** 30))

    sda = SDA(n_vis = 48**3, np_rng = np_rng)
    sda.add_DA(100)
    sda.add_DA(50)
    sda.add_DA(10)

    t_x = T.matrix('x')
    t_h = sda.t_encode(t_x)

def test_op(t_x):
    t_y = t_x + 1
    return t_y
    
def test_stack(t_x):
    h = t_x
    for i in xrange(3):
        h = test_op(h)
    return h

if __name__ == '__main__':
    test_sda()
