import os
import os.path as pt
import sys
import time
import pdb
import numpy as np
import cPickle
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
theano.config.floatX = 'float32'

#from da import DA
import da

#import hlp
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
        np_rnd = np.random.RandomState(np_rnd)
 
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

        _da = da.DA(
            np_rnd = self.np_rnd,
            th_rnd = self.th_rnd,
            n_vis = n_vis,
            n_hid = n_hid)

        _da.tag = "{0:d}:{1:d}-{2:d}".format(idx, n_vis, n_hid)
        _da.idx = len(self)
        return _da

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
        
    def f_train(self, x, y = None, corrupt = 0.2, rate = 0.1,
                  lyr = None, ec = None, dc = None):

        t_x, t_y = t_hlp.wrap_shared(x, y)

        ## request unsupervised training
        t_y = t_x if t_y is None else t_y

        x = T.matrix('x')  # a batch from t_x
        y = T.matrix('y')  # a batch from t_y

        ## corrupted input
        q = self.t_corrupt(x, corrupt)

        ## output at certian layer
        z = self.t_pipe(q, lyr, ec, dc)

        ## cross entrophy
        cost = t_hlp.cross_entrophy(y, z, axis = 1)
        
        ## squared L2 norm
        dist = t_hlp.square_l2_norm(y, z, axis = 1)

        parm = self.__get_parms__(lyr, ec, dc)

        grad = T.grad(cost, parm)

        diff = [(p, p - rate * g) for p, g in zip(parm, grad)]

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

    def f_encode(self, lyr = None, dp = None):
        x = T.matrix('x')
        z = self.t_pipe(x, lyr, dp, 0)
        return theano.function([x], z, name = "SDA_encode")

def pre_train(sda, dat, rate = 0.1, epoch = 25):
    N = dat.shape[0]
    M = dat.shape[1]

    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    ## initialize lower encoding output as raw input
    x = dat               
    for i in xrange(len(sda)):
        ## compile function for just the i_th. encoder
        ## plug in previous layer's output as input
        corr = 0.2 / (i + 1)
        train = sda.f_train(
            x = x, corrupt = 0.2, rate = rate,
            lyr = i, ec = 1, dc = 1)

        start_time = time.clock()
        # go through training epochs
        for ep in xrange(epoch):
            # go through trainng set
            cost, dist = [], []             # cost, dist
            for i_batch in xrange(n_batch):
                c, d = train(i_batch * s_batch, (i_batch + 1) * s_batch)
                cost.append(c)
                dist.append(d)
            print 'Training ep {}, cost {}, dist {}'.format(
                ep, np.mean(cost), np.mean(dist))
            sys.stdout.flush()

        end_time = time.clock()
        training_time = end_time - start_time
        print 'ran for {:.2f}m'.format(training_time / 60.)

        ## prepare input for next layer's training
        x = sda.t_encode(x, ly = i, dp = 1).eval()

def fine_tune(sda, dat, rate = 0.05, epoch = 20):
    N = dat.shape[0]
    M = dat.shape[1]

    x = dat
    # compute number of minibatches for training, validation and testing
    s_batch = 20
    n_batch = N / s_batch

    ## -------- TRAINING --------
    train = sda.f_train(x = x, corrupt = 0.05, rate = rate)

    start_time = time.clock()
    # go through training epochs
    for ep in xrange(epoch):
        # go through trainng set
        cost, dist = [], []             # cost, dist
        for i_batch in xrange(n_batch):
            c, d = train(i_batch * s_batch, (i_batch + 1) * s_batch)
            cost.append(c)
            dist.append(d)
        print 'Training ep {}, cost {}, dist {}'.format(
            ep, np.mean(cost), np.mean(dist))
        sys.stdout.flush()
        
    end_time = time.clock()
    training_time = (end_time - start_time)
    print 'ran for {:.2f}m'.format(training_time / 60.)

def save(fo, s):
    import gzip
    with gzip.open(fo, 'wb') as gz:
        cPickle.dump(s, gz, cPickle.HIGHEST_PROTOCOL)

def load(fi):
    import gzip
    with gzip.open(fi, 'rb') as gz:
        return cPickle.load(gz)
    
def build_sda(dat, seed = None):
    N = dat.shape[0]
    M = dat.shape[1]
    H = [M*2]
    while(H[-1] > 20):
        H.append(H[-1]/4)
    return SDA(n_vis = M, n_hid = H, np_rnd = seed)
    
def train_sda(sda, dat):
    ## load wm vertex sample
    N = dat.shape[0]
    M = dat.shape[1]

    print "pre-train:"
    pre_train(sda, dat, rate = 0.10, epoch = 25)

    print "fine-tune 1, rate = 0.05:"
    fine_tune(sda, dat, rate = 0.05, epoch = 30)

    # print "fine-tune 2, rate = 0.02:"
    fine_tune(sda, dat, rate = 0.02, epoch = 40)

def test():
    ## load wm vertex sample
    dat = np.load('tmp/rh11DA5.npz')['vtx'][['xyz', 'tck']]
    from itertools import chain
    xyz = [t_hlp.rescale01(dat['xyz'][a]) for a in 'xyz']
    sta = [t_hlp.rescale01(dat[a]) for a in ['tck']]
    dat = np.hstack(chain(xyz, sta))
    sda = build_sda(dat)
    train_sda(sda, dat)
    return sda, dat

def work(tsk):
    ## load data
    dat = np.load(pt.join(tsk['src'], tsk['wms'] + '.npz'))
    sbj, vtx = dat['sbj'], dat['vtx']

    ## extract xyz and thickness, and flatten these features
    x = [t_hlp.rescale01(vtx['xyz'][a]) for a in 'xyz']

    ## arrange surface thickness
    if np.count_nonzero(vtx['tck']) / vtx['tck'].size < 0.9:
        print "xt: proportio of zero exceed 0.1 in thickness"
        return
    x.append(t_hlp.rescale01(vtx['tck']))
    x = np.hstack(x)

    ## train SDA
    if not tsk.has_key('sda'):
        sda = build_sda(x, seed = tsk['seed'])
        train_sda(sda, x)
        tsk['sda'] = sda

    ## encode the features into low dimensional form
    if not tsk.has_key('enc'):
        enc = sda.f_encode()(x)
        tsk['enc'] = enc

    ## white matter surface sample id
    wms = tsk['wms']

    ## save binary: SDA, vertices, encodings and subjects
    save(pt.join(tsk['dst'], wms + '.pgz'), tsk)
    
    ## export encoding in ascii
    with open(pt.join(tsk['xpt'], wms + '.enc'), 'wb') as f:
        np.savetxt(f, enc, delimiter = '\t')

    ## export subjects in ascii
    with open(pt.join(tsk['xpt'], wms + '.sbj'), 'wb') as f:
        np.savetxt(f, sbj, fmt = "%s", delimiter = '\t')

    print "xt: success"
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as pk:
            tsk = cPickle.load(pk)
        work(tsk)
    else:
        tsk = '../../hpc/trained_sda/09_0078/script/lh05F05.pk'
        with open(tsk) as pk:
            tsk = cPickle.load(pk)
        pass
