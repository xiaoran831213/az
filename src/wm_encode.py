import os.path as pt
import pdb
import numpy as np
import rdm
import rdm.hlp as hlp
import hlp as rut_hlp
from glob import glob as gg
from rdm.sae import SAE

## fine turning
def work_3(fi, fo, dpt = [0, 4], ovr = 0):
    ## load data output from step 1: pre-training
    dat = rut_hlp.load_pgz(fi)

    nnt = dat['nnt']
    vt, wm = dat['vt'], dat['nm']

    ## save binary: SDA, vertices, encodings and subjects
    if pt.isfile(fo) and ovr <2:
        print "update:", fo; sys.stdout.flush()
    else:
        print 'create:', fo; sys.stdout.flush()

    from collections import OrderedDict as OD
    enc = OD()
    for f, d in [(n[0], d) for n in nnt for d in dpt]:
        print 'encode', wm, f, d; sys.stdout.flush()
        c0 = hlp.rescale01(vt[f])
        if d > 0:
            enc[f, d] = nnt[f, 'stk'].sub(d).ec(c0).eval()
        else:
            enc[f, 0] = c0

    dat['enc'] = enc

    ## save python data
    rut_hlp.save_pgz(fo, dat)
    print 'saved:', fo; sys.stdout.flush()
    print "xt: success"

def main(src, dst):
    for fi in gg(pt.join(src, '*.pgz')):
        fo = pt.join(dst, pt.basename(fi))
        work_3(fi, fo)
    
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
        src = pt.expandvars('$AZ_AFTN')

    if len(sys.argv) > 3:
        dst = sys.argv[3]
    else:
        dst = pt.expandvars('$AZ_AENC')

    if wms is not None:
        main(wms, src, dst)
