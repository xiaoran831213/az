import os.path as pt
import pdb
import numpy as np
import hlp
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

    print 'pre-train:', stk.dim
    x = dat
    r = rate
    for ae in stk:
        print ae.dim
        t = Trainer(ae.z, src = x, xpt = x, lrt = r)
        t.tune(epoch, 10)
        x = ae.y(x).eval()
        r = r * 2.0
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

def fine_tune(stk, dat, rate = 0.01, epoch = 50):
    print 'find-tune:', stk.dim
    x = dat
    dpt = len(stk)
    t = Trainer(stk.z, src = x, xpt = x, lrt = rate/dpt)
    t.tune(epoch * dpt, 10)
    timer = time.clock() - timer
    print 'ran for {:.2f}m\n'.format(timer / 60.)

def work(tsk, ftr = ['slc', 'tck'], eph = 100, ovr = 0):
    ## load data
    dst, src, wms = tsk['dst'], tsk['src'], tsk['wms']
    dat = np.load(pt.join(src, wms + '.npz'))
    sbj, vtx = dat['sbj'].tolist(), dat['vtx'][ftr]
    del dat

    ## save binary: SDA, vertices, encodings and subjects
    fo = pt.join(dst, wms + '.pgz')
    if pt.isfile(fo):
        if ovr < 2:
            print "exists:", fo
            for k, v in rut_hlp.load_pgz(fo).iteritems():
                tsk[k] = v
        else:
            print "overwrite:", fo

    ## quality check
    for fn, fv in [(fn, vtx[fn]) for fn in ftr]:
        if np.count_nonzero(fv) / float(fv.size) > 0.9:
            continue
        print "xt: 0s exceed 10% in {}/{}['{}']".format(
            src, wms, fn)
        return

    ## get or create network dictionary
    if not tsk.has_key('nnt'):
        tsk['nnt'] = {}
    nnt = tsk['nnt']

    if not tsk.has_key('enc'):
        tsk['enc'] = {}
    enc = tsk['enc']
    
    ## train each feature seperately for now
    ## fn: feature name, dd: power of dimension divisor
    print 'wm surface: ', wms
    from itertools import product
    for fn, dd in product(ftr, xrange(1, 1 + 10)):
        ## decide input value and dimensions
        fv = hlp.rescale01(vtx[fn])
        dm = fv.shape[-1]
        if dm / 2** dd < 8:
            continue

        ## use existing or create new network
        dm = [dm / 2**d for d in xrange(1 + dd)]

        ## re-train or non-exist
        if not nnt.has_key((fn, dd)):        # new network or re-train
            nt = SAE(dm)
            nnt[(fn, dd)] = nt
        elif ovr > 0:                                       # continue
            nt = nnt[(fn, dd)]
        else:                                          # skip existing
            print "xt: {}.{} exists.".format(fn, dd)
            continue

        ## train the network
        print 'feature: ', fn
        train(nt, fv, rate = 0.01, epoch = eph)

        ## encode the feature
        enc[(fn, dd)] = nt.ec(fv).eval()
        
    ## save
    rut_hlp.save_pgz(fo, tsk)

def rexp(tsk):
    ## export encoded surface to R
    import rpy2
    import rpy2.robjects as robjs
    import rpy2.rlike.container as rlc
    R = rpy2.robjects.r

    dst, src, wms = tsk['dst'], tsk['src'], tsk['wms']
    enc = tsk['enc']

    rds = R.readRDS(pt.join(src, wms + '.rds'))
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

    x = np.load(pt.expandvars('$AZ_SP1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2, d/4, d/8, d/16]
    m = SAE(dim=dim)
    return x, m

if __name__ == '__main__':
    import os
    import sys
    hlp.set_seed(120)

    ## add project root to python path
    if not os.environ['AZ_PRJ'] in sys.path:
        sys.path.insert(0, os.environ['AZ_PRJ'])
    import src.hlp as rut_hlp
    
    ## parse arguments
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
                work(tsk, ovr = 2, eph = 100)
                
            tsk = pt.expandvars('$AZ_EC1/tsk/rh0763B.pk')
            with open(tsk) as pk:
                tsk = cPickle.load(pk)
                work(tsk, ovr = 2, eph = 100)
