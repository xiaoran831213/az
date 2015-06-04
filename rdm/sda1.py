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
import pdb

def __get_name__(t_x):
    if hasattr(t_x, 'name'):
        name = getattr(t_x, 'name')
    else:
        name = "x:{:d}".format(t_x.size / t_x.shape[0])
    if name == None:
        name = ""
    return name

def __set_name__(t_x, name):
    if hasattr(t_x, 'name'):
        setattr(t_x, 'name', name)
        
# start-snippet-1
class SDA(list):
    """ Stacked denoising auto-encoder class (SdA) """
    def __init__(self, n_vis, n_hid = None, np_rnd = None, th_rnd = None):
        """ This class is made to support a variable number of layers. """
        if not np_rnd:
            np_rnd = np.random.RandomState(123)
 
        if not th_rnd:
            th_rnd = RandomStreams(np_rnd.randint(2 ** 30))

        self.n_vis = n_vis                # number of visible feature
        self.np_rnd = np_rnd
        self.th_rnd = th_rnd

        if not n_hid:
            n_hid = ()
        if hasattr(n_hid, '__len__'):
            self.extend(n_hid)
        else:
            self.append(n_hid)

    def __mk_da__(self, n_hid):
        idx = len(self)
        if idx:
            n_vis = self[-1].n_hid
        else:
            n_vis = self.n_vis

        da = da1.DA(
            np_rnd = self.np_rnd,
            th_rnd = self.th_rnd,
            n_vis = n_vis,
            n_hid = n_hid)
        da.tag = "{0:d}:{1:d}-{2:d}".format(idx, n_vis, n_hid)
        da.idx = len(self)
        return da

    ## override list.extend
    def extend(self, n_hds):
        return super(SDA, self).extend(
            self.__mk_da__(n_hid) for n_hid in n_hds)

    ## override list.append
    def append(self, n_hid):
        return super(SDA, self).append(
            self.__mk_da__(n_hid))
        
    def __get_stack__ (self, ly = None, ec = None, dc = None):
        if ly == None:
            ly = 0

        if ly < 0:
            ly = len(self) + ly

        if ec == None:
            lh = len(self)
        else:
            lh = min(ly + ec, len(self))

        if dc == None:
            lz = 0
        else:
            lz = max(lh - dc, 0)

        ec = self[ly:lh]
        dc = self[lz:lh]
        dc.reverse()
        return ec, dc
        
    def __get_parms__(self, ly = 0, ec = None, dc = None):
        ec, dc = self.__get_stack__(ly, ec, dc)
        parms = []
        for da in ec:
            parms.append(da.t_w)
            parms.append(da.t_b)
        for da in dc:
            parms.append(da.t_b_prime)
        return parms

    def t_pipe(self, t_x, ly = 0, ec = None, dc = None):
        ec, dc = self.__get_stack__(ly, ec, dc)
        name = __get_name__(t_x)
        
        ## build pipe expression
        for da in ec:
            t_x = da.t_encode(t_x)
            name += "|{0:d}:{1:d}-{2:d}".format(da.idx, da.n_vis, da.n_hid)
        for da in dc:
            t_x = da.t_decode(t_x)
            name += "|{0:d}:{2:d}-{1:d}".format(da.idx, da.n_vis, da.n_hid)

        __set_name__(t_x, name)
        return t_x

    def t_encode(self, t_x, ly = 0, dp = None):
        return self.t_pipe(t_x, ly, ec = dp, dc = 0)

    def t_decode(self, t_x, ly = None, dp = None):
        if ly == None:
            ly = len(self)
        return self.t_pipe(t_x, ly, ec = 0, dc = dp)

    def t_corrupt(self, t_x, lvl):
        return self.th_rnd.binomial(
            size = t_x.shape, n = 1, p = 1 - lvl,
            dtype = T.config.floatX) * t_x

    def f_train(self, t_x, t_y = None, t_corrupt = 0.2, t_rate = 0.1,
                  ly = None, ec = None, dc = None):

        if t_y == None:    # unsupervised training
            t_y = t_x
        x = T.matrix('x')  # a batch from t_x
        y = T.matrix('y')  # a batch from t_y
        q = self.t_corrupt(x, t_corrupt)
        z = self.t_pipe(q, ly, ec, dc)

        L = - T.sum(y * T.log(z) + (1 - y) * T.log(1 - z), axis=1)
        cost = T.mean(L)

        dist = T.mean(T.sqrt(T.sum((x - z) ** 2, axis = 1)))

        parm = self.__get_parms__(ly, ec, dc)

        grad = T.grad(cost, parm)

        diff = [(p, p - t_rate * g) for p, g in zip(parm, grad)]

        t_fr = T.iscalar()
        t_to = T.iscalar()
        return theano.function(
            [t_fr, t_to],
            [cost, dist],
            updates = diff,
            givens = {x : t_x[t_fr:t_to], y : t_y[t_fr:t_to]},
            name = "SDA_trainer")

    def f_pred(self, ly = None, ec = None, dc = None):
        x = T.matrix('x')
        z = self.t_pipe(x, ly, ec, dc)
        return theano.function([x], z, name = "SDA_pred")

def test_sda(dat):
    x = hlp.get_pk(dat)
    x = x.reshape(x.shape[0], -1)
    sda = SDA(n_vis = x.size/x.shape[0], n_hid = (270, 90, 30, 10))

    print "pre-train:"
    test_pre_train(sda, x)
    z = sda.f_pred()(x)
    print "AUC after pre-train: ", hlp.AUC(x, z)

    print "fine-tune:"
    test_fine_tune(sda, x)
    z = sda.f_pred()(x)
    print "AUC after fine-tune: ", hlp.AUC(x, z)

    print hlp.AUC(x, z)
    return sda
    
def test_pre_train(sda, dat):
    N = dat.shape[0]
    M = dat.size/dat.shape[1]
    s_x = theano.shared(np.asarray(
        dat, dtype = theano.config.floatX), borrow = True)

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    for i in xrange(len(sda)):
        t_i = sda.t_encode(s_x, ly = 0, dp = i)
        corrupt = 0.2
        train = sda.f_train(
            t_x = t_i, t_corrupt = 0.2, t_rate = 0.1,
            ly = i, ec = 1, dc = 1)

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

def test_fine_tune(sda, dat):
    N = dat.shape[0]
    M = dat.size/dat.shape[1]
    s_x = theano.shared(np.asarray(
        dat, dtype = theano.config.floatX), borrow = True)

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    train = sda.f_train(
        t_x = s_x, t_corrupt = 0.2, t_rate = 0.1)

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
    da = da1.DA(np_rnd = np_rng, n_vis = M, n_hid = 90)

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
