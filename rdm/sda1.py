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

import hlp

# start-snippet-1
class SDA(list):
    """ Stacked denoising auto-encoder class (SdA) """

    def __init__(self, n_vis, np_rng, th_rng=None):
        """ This class is made to support a variable number of layers. """
        if not th_rng:
            th_rng = RandomStreams(np_rng.randint(2 ** 30))

        self.n_vis = n_vis                # number of visible feature
        self.np_rng = np_rng
        self.th_rng = th_rng
        

    def append(self, n_hid):
        if len(self):
            n_vis = self[-1].n_hid
        else:
            n_vis = self.n_vis

        da = da1.DA(
            np_rng = self.np_rng,
            th_rng = self.th_rng,
            n_vis = n_vis,
            n_hid = n_hid)
        return super(SDA, self).append(da)
        
    def __get_stack__ (self, lyr = None, dp_enc = None, dp_dec = None):
        if lyr == None:
            lyr = 0
            dp_enc = -1
            dp_dec = -1
        if lyr < 0:
            lyr = len(self) - lyr
        if dp_enc == None:
            dp_enc = 1
        if dp_enc < 0:
            dp_enc = sys.maxint

        if dp_dec == None:
            dp_dec = 1
        if dp_dec < 0:
            dp_dec = sys.maxint

        l_h = min(lyr + dp_enc, len(self))  # last hidden layer + 1
        l_z = max(l_h - dp_dec, 0)               # last re-con layer
        ec = self[lyr:l_h]
        dc = self[l_z:l_h]
        dc.reverse()
        return ec, dc
        
    def __get_parms__(self, lyr = None, dp_enc = None, dp_dec = None):
        ec, dc = self.__get_stack__(lyr, dp_enc, dp_dec)
        parms = []
        for da in ec:
            parms.append(da.t_w)
            parms.append(da.t_b)
        for da in dc:
            parms.append(da.t_b_prime)
        return parms

    def t_pipe(self, t_x, lyr = None, dp_enc = None, dp_dec = None):
        ec, dc = self.__get_stack__(lyr, dp_enc, dp_dec)
        for da in ec:
            t_x = da.t_encode(t_x)
        for da in dc:
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

    def f_train(self, t_data, t_corrupt = 0.2, t_rate = 0.1,
                  lyr = None, dp_enc = None, dp_dec = None):
        x = T.matrix('x')
        q = self.t_corrupt(x, t_corrupt)
        z = self.t_pipe(q, lyr, dp_enc, dp_dec)

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
            name = "SDA_trainer")

    def f_pred(self, lyr = None, dp_enc = None, dp_dec = None):
        x = T.matrix('x')
        z = self.t_pipe(x, lyr, dp_enc, dp_dec)
        return theano.function([x], z, name = "SDA_pred")

def test_sda():
    dat = hlp.get_pk('dat/d48/1003')
    N = dat.shape[0]
    dat = np.reshape(dat, (N, -1))
    M = dat.shape[1]
    s_x = theano.shared(np.asarray(
        dat, dtype = theano.config.floatX), borrow = True)

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    np_rng = np.random.RandomState(123)
    sda = SDA(n_vis = M, np_rng = np_rng)
    sda.append(100)
    sda.append(100)
    sda.append(50)

    ## -------- TRAINING --------
    train = sda.f_train(
        t_data = s_x, t_corrupt = 0.2, t_rate = 0.1,
        lyr = 0, dp_enc = 3, dp_dec = 3)

    start_time = time.clock()
    # go through training epochs
    for epoch in xrange(15):
        # go through trainng set
        c, d = [], []                     # cost, dist
        for i_batch in xrange(n_batch):
            r = train(i_batch * s_batch, (i_batch + 1) * s_batch)
            c.append(r[0])
            d.append(r[1])
        print 'Training epoch %d, cost %f, dist %f' % (epoch, np.mean(c), np.mean(d))

    end_time = time.clock()
    training_time = (end_time - start_time)
    print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))
    return sda

def test_da1():
    dat = hlp.get_pk('dat/d48/1003')
    N = dat.shape[0]
    dat = np.reshape(dat, (N, -1))
    M = dat.shape[1]
    s_x = theano.shared(np.asarray(
        dat, dtype = theano.config.floatX), borrow = True)

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    np_rng = np.random.RandomState(123)
    da = da1.DA(np_rng = np_rng, n_vis = M, n_hid = 100)

    # ## -------- TRAINING --------
    train = da.f_train(t_x = s_x, t_corrupt = 0.2, t_rate = 0.1)
    ## we know S_x.eval().shape[0] = 48**3
    start_time = time.clock()
    # go through training epochs
    for epoch in xrange(3):
        # go through trainng set
        c, d = [], []                     # cost, dist
        for i_batch in xrange(n_batch):
            r = train(i_batch * s_batch, (i_batch + 1) * s_batch)
            c.append(r[0])
            d.append(r[1])
        print 'Training epoch %d, cost %f, dist %f' % (epoch, np.mean(c), np.mean(d))

    end_time = time.clock()
    training_time = (end_time - start_time)
    print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))
    return da
    
if __name__ == '__main__':
    pass
