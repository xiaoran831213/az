import os
import os.path as pt
import pdb
import numpy as np
import hlp
from hlp import T
import trainer
from trainer import Trainer
import sae
from sae import SAE

def train(stk, dat, rate = 0.01, epoch = 50):
    """
    pre-train each auto encoder in the stack
    """
    import time
    timer = time.clock()

    print 'pre-train:'
    x = dat
    r = rate
    for ae in stk:
        t = Trainer(ae.z, src = x, xpt = x, lrt = r)
        t.tune(epoch, 10)
        x = ae.y(x).eval()
        r = r * 2.0
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

    print 'find-tune:'
    x = dat
    dpt = len(stk)
    t = Trainer(stk.z, src = x, xpt = x, lrt = rate/dpt)
    t.tune(epoch * dpt, 10)
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

def work(tsk):
    ## load data
    dst, src, wms = tsk['dst'], tsk['src'], tsk['wms']
    dat = np.load(pt.join(src, wms + '.npz'))

    ##  pick xyz and thickness for now
    sbj, vtx = dat['sbj'].tolist(), dat['vtx']
    del dat

    ## save binary: SDA, vertices, encodings and subjects
    fo = pt.join(dst, wms + '.pgz')
    if pt.isfile(fo):
        print fo + ": exists"
        tsk = hlp.load_pgz(fo)

    ## quality check
    for fn, fv in [(fn, vtx[fn]) for fn in vtx.dtype.names]:
        if np.count_nonzero(fv) / float(fv.size) > 0.9:
            continue
        print "xt: 0s exceed 10% in {}/{}['{}']".format(
            src, wms, fn)
        return
    
    ## the dimension divisor
    div = ([2**j for j in xrange(1, i)] for i in xrange(2, 12))

    from itertools import product
    from math import log

    ## get or create network dictionary
    if not tsk.has_key('nnt'):
        tsk['nnt'] = {}
    nnt = tsk['nnt']

    if not tsk.has_key('enc'):
        tsk['enc'] = {}
    enc = tsk['enc']
    
    ## train each feature seperately for now
    ## fn: feature name, dd: power of dimension divisor
    for fn, dd in product(vtx.dtype.names, xrange(1, 1 + 10)):
        ## decide input value and dimensions
        fv = hlp.rescale01(vtx[fn])
        dm = fv.shape[-1]
        if dm / 2** dd < 4 :
            continue
        
        dm = [dm / 2**d for d in xrange(1 + dd)]
        
        ## get or create the neural network
        if not nnt.has_key((fn, dd)):
            nnt[(fn, dd)] = SAE(dm)
        nt = nnt[(fn, dd)]

        ## train the network
        train(nt, fv, rate = 0.01, epoch = 100)

        ## encode the feature
        enc[(fn, dd)] = nt.ec(fv).eval()
        
    ## save
    hlp.save_pgz(fo, tsk)
    del nnt, tsk

    # ## export data to R
    # import rpy2
    # import rpy2.robjects as robjs
    # R = rpy2.robjects.r

    # fo = pt.join(dst, wms + '.rds')
    # if pt.isfile(fo):
    #     print fo + ": exists"
    # else:
    #     ## 1) raw input, (f -- feature, v -- vertex)
    #     x = [f for v in vtx.flat for f in v.item()]
    #     x = R.array(x, dim = [len(ftr), vtx.shape[1], vtx.shape[0]], dimnames = [ftr, [], sbj])
    #     x = R.aperm(x, [2,1,3])
        
    #     ## 2) encoding, c: code, s: subject
    #     y = enc.flatten().tolist()
    #     y = R.array(y, dim = [enc.shape[1], enc.shape[0]], dimnames = [[], sbj])
    #     y = R.aperm(y, [2,1])

    #     ## group R data in a environment and save to binary
    #     e = robjs.Environment()
    #     e['x'] = x
    #     e['y'] = y
    #     R.saveRDS(e, fo)

    print "xt: success"
    
def test_sae():
    hlp.set_seed(120)

    x = np.load(pt.expandvars('$AZ_SP1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2, d/4, d/8, d/16]
    m = SAE(dim=dim)
    return x, m
    
if __name__ == '__main__':
    import sys
    import cPickle
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as pk:
            tsk = cPickle.load(pk)
        work(tsk)
    else:
        if not 'tsk' in dir():
            tsk = pt.expandvars('$AZ_EC1/tsk/lh001F1.pk')
            with open(tsk) as pk:
                tsk = cPickle.load(pk)
