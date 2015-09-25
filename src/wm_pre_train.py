import os.path as pt
import pdb
import numpy as np
import rdm
import rdm.hlp as hlp
import hlp as rut_hlp
from rdm.sae import SAE
from rdm.trainer import Trainer

def pre_train(stk, dat, rate = 0.01, epoch = 100):
    """
    pre-train each auto encoder in the stack
    """
    import time
    tm = time.clock()

    print 'pre-train:', stk.dim; sys.stdout.flush()
    x = dat
    r = rate
    d0 = dat.shape[1]
    for ae in stk.sa:
        print ae.dim
        t = Trainer(ae, src = x, dst = x, lrt = r)
        t.tune(epoch, 10)
        x = ae.ec(x).eval()
        ## adjust movement speed by parameter count
        r = rate * d0 / x.shape[1]
        
    tm = time.clock() - tm
    print 'ran for {:.2f}m\n'.format(tm/60.); sys.stdout.flush()

def work_1(fi, fo, ftr = 'tck', ep = 100, lr = 0.01, ovr = 0):
    ## load data
    dat = np.load(fi)
    vtx, wms = hlp.rescale01(dat['vt'][ftr]), dat['nm']

    ## dimensions
    dim = [dat['vi'].shape[0] / (2**i) for i in xrange(8)]

    ## load binary: SDA, vertices, subjects
    if pt.isfile(fo) and ovr < 2:
        print "update:", fo; sys.stdout.flush()
        tsk = rut_hlp.load_pgz(fo)
    else:
        print 'create:', fo
        tsk = dict([(k, v) for k, v in dat.iteritems()])
        tsk['nnt'], tsk['eph'] = {}, {}
        tsk['dim'] = dim
    ## close npz file pointer.
    dat.close()

    ## get network and encode dictionary
    nnt, eph = tsk['nnt'], tsk['eph']

    ## pre train
    print 'surface: ', wms; sys.stdout.flush()

    print 'feature:', ftr
    key = (ftr, 'stk')
    if not nnt.has_key(key):
        nnt[key] = SAE.from_dim(dim)
        eph[key] = 0
        print "fill: {}.{}".format(ftr, 'stk')
        ovr = 2
    elif ovr == 0:
        print "skip: {}.{}".format(ftr, 'stk')
    else:
        print "more: {}.{}".format(ftr, 'stk')

    if ovr > 0:
        pre_train(nnt[key], vtx, rate = lr, epoch = ep)
        eph[key] += ep
    
    ## save python data
    rut_hlp.save_pgz(fo, tsk)
    print 'pre-trained:', fo; sys.stdout.flush()

def main(wms, src, dst):
    fi = pt.join(src, wms + '.npz')
    fo = pt.join(dst, wms + '.pgz')
    work_1(fi, fo, ep = 100, lr = 0.01, ovr = 1)

if __name__ == '__main__':
    import os
    import sys
    hlp.set_seed(None)

    if len(sys.argv) > 1:
        wms = sys.argv[1]
    else:
        wms = None
        
    if len(sys.argv) > 2:
        src = sys.argv[2]
    else:
        src = pt.expandvars('$AZ_APEC')

    if len(sys.argv) > 3:
        dst = sys.argv[3]
    else:
        dst = pt.expandvars('$AZ_APTN')

    if wms is not None:
        main(wms, src, dst)
