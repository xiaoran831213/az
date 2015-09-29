import os.path as pt
import pdb
import numpy as np
import rdm
import rdm.hlp as hlp
import hlp as rut_hlp
from rdm.sae import SAE
from rdm.trainer import Trainer

## fine turning
def work_2(fi, fo, ftr = 'tck', dpt = None, ep = 100, lr = 0.01, ovr = 0):
    ## load data output from step 1: pre-training
    dat = rut_hlp.load_pgz(fi)
    vt, wm = hlp.rescale01(dat['vt'][ftr]), dat['nm']
    key = (ftr, 'stk')

    ## save binary: SDA, vertices, encodings and subjects
    if pt.isfile(fo) and ovr <2:
        print "update:", fo; sys.stdout.flush()
        tsk = rut_hlp.load_pgz(fo)
        if not tsk['nnt'].has_key(key):
            tsk['nnt'][key] = dat['nnt'][key]
    else:
        tsk = dat
        print 'create:', fo

    ## get network and encode dictionary
    stk, eph = tsk['nnt'][key], tsk['eph']

    if dpt is not None:
        stk = stk.sub(dpt)
        
    ## fine tune
    print 'surface & feature: ', wm, ftr; sys.stdout.flush()

    ## learning rate is adjusted by # of parameters
    lr = lr / sum(tsk['dim'])

    import time
    tm = time.clock()
    
    print 'find-tune:', stk.dim; sys.stdout.flush()
    t = Trainer(stk, src = vt, dst = vt, lrt = lr)
    t.tune(ep, 10)
    tm = time.clock() - tm
    print 'ran for {:.2f}m\n'.format(tm/60.);  sys.stdout.flush()

    if(t.ecst() < 0):
        print "xt: failure"
    else:
        ## save python data
        rut_hlp.save_pgz(fo, tsk)
        print 'saved:', fo; sys.stdout.flush()
        print "xt: success"

def main(wms, src, dst):
    fi = pt.join(src, wms + '.pgz')
    fo = pt.join(dst, wms + '.pgz')
    work_2(fi, fo, dpt = 4, ep = 40, lr = 0.005, ovr = 1)

if __name__ == '__main__':
    import os
    import sys
    if len(sys.argv) > 1:
        wms = sys.argv[1]
    else:
        wms = None
        
    if len(sys.argv) > 2:
        src = sys.argv[2]
    else:
        src = pt.expandvars('$AZ_APTN')

    if len(sys.argv) > 3:
        dst = sys.argv[3]
    else:
        dst = pt.expandvars('$AZ_AFTN')

    if wms is not None:
        main(wms, src, dst)
