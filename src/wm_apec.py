#!/usr/bin/env python
import pdb
import numpy as np
import sys
import os
import os.path as pt
from glob import glob as gg

def apec(src, dst, sn, hm, ovr=0):
    """
    given a directory of subjects, extract anatomical regions
    src: subject source
    dst: where to put extracted regions
    sn: region serial number, 2~35
    hm: hemesphere, lh or rh
    """
    ## output path
    fo = pt.join(dst, '{}{:02d}.npz'.format(hm, sn))
    if pt.exists(fo) and not ovr:
        print 'exists:', fo
        return
    
    ## fetch anatomical peceration table
    at = np.load('apec.npz')

    ## vertex indices
    vi = (at[hm]==sn).nonzero()[0]

    ## anatomy region id and name
    id, nm = at['tb'][sn]['id'], at['tb'][sn]['nm']

    ## surface, and subject index
    vt, sb = [], []

    ## iterate all subjects
    print 'xt: extract', hm, at['tb'][sn]['nm'], 'from ', src, ':'

    for fn in gg(pt.join(src, '*.npz')):
        sb.append(pt.basename(fn).split('.')[0])
        print sb[-1]
        sys.stdout.flush()
        wm = np.load(fn)

        ## extract surfaces for subject {sb}
        ## hm: hemisphere, vi: vertex indices
        vt.append(wm[hm][vi])
        wm.close()

    ## write the samples to file in numpy format.
    print 'xt: write surface to ', dst
    sys.stdout.flush()

    vt=np.vstack(vt)
    sb=np.array(sb)
    np.savez_compressed(
        fo, sb=sb, vt=vt, vi=vi, hm=hm, sn=sn, id=id, nm=nm)
    
    print 'xt: success'
    sys.stdout.flush()

def main():
    import sys
    src = pt.expandvars('$AZ_AVTX')
    dst = pt.expandvars('$AZ_APEC')
    for sn in xrange(3, 37):
        apec(src, dst, sn, 'lh')
        apec(src, dst, sn, 'rh')

if __name__ == "__main__":
    main()
