import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import cPickle
from itertools import izip

def __wm_neighbor__(nb, cv, sz):
    """
    According to vertex neighborhood {nb}, sample {sz}
    neighbors around the center vertex {cv}
    return vertex indices including the center one
    """
    ## mark vertices
    idx = [cv]
    mrk = set(idx)
    ## neighbor format: v_idx, n_nbr, nb[0], nb[1], ... nb[n_nbr-1]
    for i in xrange(sz):
        if len(idx) < sz:
            new = set(nb[idx[i]][2:]).difference(mrk)
            mrk.update(new)
            idx.extend(new)
        else:
            return idx[:sz]

def __sample_wm__(wrk):
    """
    given a list of center vertices {cvs}, and hemispheres {hms},
    pick regions from WM surface across subject in {src}
    """
    ## fetch working specifications
    src = wrk['src']     # source directory with subjects
    dst = wrk['dst']     # target filenames
    hms = wrk['hms']     # hemispheres
    cvs = wrk['cvs']     # center vertices
    nbs = wrk['nbs']     # neighborhoood
    sz = wrk['sz']       # region size

    ## gather lists of neighbor indices
    lfo, lvi = [], []
    for hm, nb, cv in izip(hms, nbs, cvs):
        lvi.append(__wm_neighbor__(nb, cv, sz))
        lfo.append(pt.join(dst, '{}{:05X}.npy'.format(hm, cv)))
        
    ## make lists of sampled surfaces and subject IDs
    lsf = [[]] * len(lvi)
    print 'xt: sample ', len(lvi), 'WM areas from ', src, ':'
    for fn in gg(pt.join(src, '*')):
        if not fn.endswith('wm.npz'):
            continue
        print pt.basename(fn)
        wm = np.load(fn)
        for si, hm, vi in izip(xrange(len(lvi)), hms, lvi):
            lsf[si].append(wm[hm][vi])

    print 'xt: write WM area samples to ', dst, ':'
    for sf, fo in izip(lsf, lfo):
        np.save(fo, np.array(sf))
        print fo, ": created"
    print 'xt: success'

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        wrk = '../tmp/wm_samples/0B_0078/script/000_00.pk'
    else:
        wrk = sys.argv[1]
    with open(wrk) as pk:
        wrk = cPickle.load(pk)
    __sample_wm__(wrk)
