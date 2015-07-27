import os
import os.path as pt
import sys
import time
import pdb
import numpy as np
import cPickle
import hlp
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
from sda import SDA
import rpy2
import rpy2.robjects as robjs
R = rpy2.robjects.r

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
    pre_train(sda, dat, rate = 0.10, epoch = 20)

    print "fine-tune 1, rate = 0.05:"
    fine_tune(sda, dat, rate = 0.05, epoch = 30)

    print "fine-tune 2, rate = 0.02:"
    fine_tune(sda, dat, rate = 0.02, epoch = 30)


def work(tsk):
    ## load data
    dst = tsk['dst']
    src = tsk['src']
    wms = tsk['wms']
    dat = np.load(pt.join(src, wms + '.npz'))

    ## only pick xyz and thickness for now
    ftr = ['x', 'y', 'z', 'tck']
    sbj, vtx = dat['sbj'].tolist(), dat['vtx'][ftr]
    del dat

    ## quality check
    if np.count_nonzero(vtx['tck']) / vtx['tck'].size < 0.9:
        print "xt: proportio of zero exceed 0.1 in thickness"
        return
    pdb.set_trace()
    ## save binary: SDA, vertices, encodings and subjects
    fo = pt.join(dst, wms + '.pgz')
    if pt.isfile(fo):
        print fo + ": exists"
        tsk = load(fo)
        enc = tsk['enc']
    else:
        ## 1) rescale the features to [0,1] and flatten them
        x = np.hstack([hlp.rescale01(vtx[a]) for a in ftr])

        ## 1) train SDA
        sda = build_sda(x, seed = tsk['seed'])
        train_sda(sda, x)
        tsk['sda'] = sda

        ## 2) encode
        enc = sda.f_encode()(x)
        tsk['enc'] = enc
        del sda, x

        ## 3) save
        save(fo, tsk)
    del tsk

    ## export data to R
    fo = pt.join(dst, wms + '.rds')
    if pt.isfile(fo):
        print fo + ": exists"
    else:
        ## 1) raw input, (f -- feature, v -- vertex)
        x = [f for v in vtx.flat for f in v.item()]
        x = R.array(x, dim = [len(ftr), vtx.shape[1], vtx.shape[0]], dimnames = [ftr, [], sbj])
        x = R.aperm(x, [2,1,3])
        
        ## 2) encoding, c: code, s: subject
        y = enc.flatten().tolist()
        y = R.array(y, dim = [enc.shape[1], enc.shape[0]], dimnames = [[], sbj])
        y = R.aperm(y, [2,1])

        ## group R data in a environment and save to binary
        e = robjs.Environment()
        e['x'] = x
        e['y'] = y
        R.saveRDS(e, fo)

    print "xt: success"
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as pk:
            tsk = cPickle.load(pk)
        work(tsk)
    else:
        if not 'tsk' in dir():
            with open('../../hpc/trained_sda/09_0078/script/lh05F05.pk') as pk:
                tsk = cPickle.load(pk)
