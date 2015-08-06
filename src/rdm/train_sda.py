import os
import os.path as pt
import sys
import time
import pdb
import numpy as np
import hlp
from hlp import T
import trainer
from trainer import Trainer
import sae
from sae import SAE
import cPickle

def train(stk, dat, rate = 0.01, epoch = 50):
    """
    pre-train each auto encoder in the stack
    """
    timer = time.clock()

    print 'pre-train:'
    x = dat
    r = rate
    for ae in stk:
        t = Trainer(ae.z, src = x, xpt = x, lrt = r)
        t.tune(epoch)
        x = ae.y(x).eval()
        r = r * 2.0
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

    print 'find-tune:'
    x = dat
    dpt = len(stk)
    t = Trainer(stk.z, src = x, xpt = x, lrt = rate/dpt)
    t.tune(epoch * dpt)
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

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
    
def work(tsk):
    ## load data
    dst = tsk['dst']
    src = tsk['src']
    wms = tsk['wms']
    dat = np.load(pt.join(src, wms + '.npz'))

    ##  pick xyz and thickness for now
    ftr = ['x', 'y', 'z', 'tck']
    sbj, vtx = dat['sbj'].tolist(), dat['vtx'][ftr]
    del dat

    ## quality check
    if np.count_nonzero(vtx['tck']) / vtx['tck'].size < 0.9:
        print "xt: proportio of zero exceed 0.1 in thickness"
        return
    
    ## save binary: SDA, vertices, encodings and subjects
    fo = pt.join(dst, wms + '.pgz')
    if pt.isfile(fo):
        print fo + ": exists"
        tsk = load(fo)

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
    import rpy2
    import rpy2.robjects as robjs
    R = rpy2.robjects.r

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
    
def test_sae():
    hlp.set_seed(120)

    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2, d/4, d/8, d/16]
    m = SAE(dim=dim)
    return x, m
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as pk:
            tsk = cPickle.load(pk)
        work(tsk)
    else:
        pass
        # if not 'tsk' in dir():
        #     with open('../../hpc/trained_sda/09_0078/script/lh05F05.pk') as pk:
        #         tsk = cPickle.load(pk)
