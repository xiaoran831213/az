import os
import sys
import time
import pdb
import numpy as np

import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams

#from da import DA
from da1 import DA
FT = np.dtype('<f4')
theano.config.floatX = FT

#import hlp
import pdb
import t_hlp

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

        da = DA(
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

    def f_encode(self, ly = 0, dp = None):
        x = T.matrix('x')
        y = self.t_encode(x, ly, dp)
        return theano.function([x], y, name = "SDA_encode")
        
    def f_train(self, x, y = None, t_corrupt = 0.2, t_rate = 0.1,
                  ly = None, ec = None, dc = None):

        t_x, t_y = t_hlp.wrap_shared(x, y)

        ## request unsupervised training
        t_y = t_x if t_y is None else t_y

        x = T.matrix('x')  # a batch from t_x
        y = T.matrix('y')  # a batch from t_y

        ## corrupted input
        q = self.t_corrupt(x, t_corrupt)

        ## output at certian layer
        z = self.t_pipe(q, ly, ec, dc)

        ## cross entrophy
        cost = t_hlp.cross_entrophy(y, z, axis = 1)
        
        ## squared L2 norm
        dist = t_hlp.square_l2_norm(y, z, axis = 1)

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

def train_sda(sda = None, x = None):
    ## load wm vertex sample
    if x is None:
        x = np.load('tmp/rh11DA5.npz')['vtx']['tck']
    if sda is None:
        sda = SDA(n_vis = x.shape[1], n_hid = (1024, 256, 64, 32, 8))

    x = np.load(dat)['vtx'][['xyz', 'tck']]
    print "pre-train:"
    pre_train(sda, x)
    #z = sda.f_pred()(x)

    # print "fine-tune:"
    # fine_tune(sda, x)
    # z = sda.f_pred()(x)

    # print hlp.AUC(x, z)
    return sda

def pre_train(sda = None, dat = None, rate = 0.1):
    ## load wm vertex sample
    if dat is None:
        dat = np.load('tmp/rh11DA5.npz')['vtx']

    tck = dat['tck']
    tck = tck.reshape(tck.shape[0], -1)
    N = dat.shape[0]
    M = dat.shape[1]
    tck = (tck - tck.min()) / (tck.max() - tck.min())

    if sda is None:
        sda = SDA(n_vis = M, n_hid = (2048, 512, 128, 32, 8))

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    x = tck               # initialize lower encoding output as raw input
    for i in xrange(len(sda)):
        ## compile function for just the i_th. encoder
        ## plug in previous layer's output as input
        train = sda.f_train(
            x = x, t_corrupt = 0.2, t_rate = rate,
            ly = i, ec = 1, dc = 1)

        start_time = time.clock()
        # go through training epochs
        for epoch in xrange(15):
            # go through trainng set
            c, d = [], []             # cost, dist
            for i_batch in xrange(n_batch):
                r = train(i_batch * s_batch, (i_batch + 1) * s_batch)
                c.append(r[0])
                d.append(r[1])
            print 'Training epoch %d, cost %f, dist %f' % (epoch, np.mean(c), np.mean(d))

        end_time = time.clock()
        training_time = (end_time - start_time)
        print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))

        ## prepare input for next layer's training
        x = sda.t_encode(x, ly = i, dp = 1).eval()

    return sda

def fine_tune(sda, dat, rate = 0.05):
    ## load wm vertex sample
    if dat is None:
        dat = np.load('tmp/rh11DA5.npz')['vtx']['tck']


    x = dat
    x = x.reshape(x.shape[0], -1)
    N = dat.shape[0]
    M = dat.shape[1]
    x = (x - x.min()) / (x.max() - x.min())

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    train = sda.f_train(
        t_x = x, t_corrupt = 0.2, t_rate = rate)

    start_time = time.clock()
    # go through training epochs
    for epoch in xrange(30):
        # go through trainng set
        c = []                     # cost, dist
        for i_batch in xrange(n_batch):
            r = train(i_batch * s_batch, (i_batch + 1) * s_batch)
            c.append(r[0])
        print 'Training epoch %d, cost %f' % (epoch, np.mean(c))
        
    end_time = time.clock()
    training_time = (end_time - start_time)
    print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))

if __name__ == '__main__':
    test_file = 'tmp/rh11DA5.npz'
    pass
