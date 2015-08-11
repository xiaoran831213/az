import os.path as pt
import pdb
import numpy as np
import hlp
import sae
from sae import SAE

def train(stk, dat, rate = 0.01, epoch = 50):
    """
    pre-train each auto encoder in the stack
    """
    import time
    from trainer import Trainer
    tm = time.clock()

    print 'pre-train:', stk.dim; sys.stdout.flush()
    x = dat
    r = rate
    for ae in stk:
        print ae.dim
        t = Trainer(ae.z, src = x, dst = x, lrt = r)
        t.tune(epoch, 10)
        x = ae.y(x).eval()
        r = r * 2.0
    tm = time.clock() - tm
    print 'ran for {:.2f}m\n'.format(tm/60.); sys.stdout.flush()

    print 'find-tune:', stk.dim; sys.stdout.flush()
    x = dat
    dpt = len(stk)

    ## the training should be slower when parameters is more numerous
    t = Trainer(stk.z, src = x, dst = x, lrt = rate/dpt)

    ## fine tune requires more steps when network goes deeper
    t.tune(epoch * dpt * 2, 10)           
    tm = time.clock() - tm
    print 'ran for {:.2f}m\n'.format(tm/60.);  sys.stdout.flush()

def __np2r__(x):
    """ Numpy array to R array """
    import rpy2
    import rpy2.robjects as robjs
    import rpy2.rlike.container as rlc

    ## R only accept list based arguments
    shp = list(x.shape)
    typ = x.dtype.name

    ## transpose is necessary because R grows the left most dimension
    ## fastest
    x = x.transpose().flatten().tolist()

    ## directly involking R.array without vector wrapping results in
    ## an array of lists of size 1
    if typ.startswith('int'):
        x = robjs.IntVector(x)
    else:
        x = robjs.FloatVector(x)

    ## now x is an array of numbers
    x = robjs.r.array(x, dim = shp)
    return x
    
def __enc2rds__(tsk):
    """ append encoded surface to R data """
    import rpy2
    import rpy2.robjects as robjs
    import rpy2.rlike.container as rlc
    R = rpy2.robjects.r

    dst, src, wms = tsk['dst'], tsk['src'], tsk['wms']

    assert(tsk.has_key('enc'))
    enc = tsk['enc']
    
    ## read source RDS
    ## 0.sfs: surfaces
    ## 1.cmx: connection matrix
    fi = pt.join(src, wms + '.rds')
    rbj = R.readRDS(fi)
    rbj.names=['sfs', 'cmx']

    ## append surface encoding
    ## 2.enc: the encoding
    from collections import OrderedDict
    od = OrderedDict((
        '{}.{}'.format(*k),
        __np2r__(enc[k]))
        for k in sorted(enc.keys()))

    enc = robjs.ListVector(od)
    rbj.rx2['enc'] = enc
    
    fo = pt.join(dst, wms + '.rds')
    R.saveRDS(rbj, fo)
    print 'saved:', fo; sys.stdout.flush()
    
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
            print "exist:", fo; sys.stdout.flush()
            for k, v in rut_hlp.load_pgz(fo).iteritems():
                tsk[k] = v
        else:
            print "renew:", fo; sys.stdout.flush()

    ## quality check
    for fn, fv in [(fn, vtx[fn]) for fn in ftr]:
        if np.count_nonzero(fv) / float(fv.size) > 0.9:
            continue
        print "xt: 0s exceed 10% in {}/{}['{}']".format(
            src, wms, fn); sys.stdout.flush()
        return

    ## get or create network and encode dictionary
    if not tsk.has_key('nnt'):
        tsk['nnt'] = {}
    nnt = tsk['nnt']

    if not tsk.has_key('enc'):
        tsk['enc'] = {}
    enc = tsk['enc']
    
    ## train each feature seperately for now
    ## fn: feature name, dd: power of dimension divisor
    print 'wm surface: ', wms; sys.stdout.flush()
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
            print "{}.{} exists.".format(fn, dd); sys.stdout.flush()
            continue

        ## train the network
        print 'feature: ', fn
        train(nt, fv, rate = 0.01, epoch = eph)

        ## encode the feature
        enc[(fn, dd)] = nt.ec(fv).eval()
        
    ## save python data
    rut_hlp.save_pgz(fo, tsk)
    print 'saved:', fo; sys.stdout.flush()

    ## append encoding to R data and save
    __enc2rds__(tsk)
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
        tsk = pt.expandvars(sys.argv[1])
        with open(tsk) as pk:
            tsk = cPickle.load(pk)
        work(tsk, eph = 20, ovr = 0)
    elif 'tsk' in dir():
        pass
    else:
        ## tsk = pt.expandvars('$AZ_EC1/tsk/lh001F1.pk')
        tsk = pt.expandvars('$AZ_EC1/tsk/rh0763B.pk')
        with open(tsk) as pk:
            tsk = cPickle.load(pk)
