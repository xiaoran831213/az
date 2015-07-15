import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import cPickle
from itertools import izip
import rpy2
import rpy2.robjects as robjs
R = rpy2.robjects.r

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

def __sample2rds__(dat, sbj, vtx, fo):
    """
    save surface samples to R binary.
    dat: samples of a surface region.
    sbj: subject list, indices the 1st dimension of {dat}
    """
    hdr = list(dat[0].dtype.names)
    ## flatten the data
    dat = [f for s in dat for v in s for f in v.item()]

    ## create 3D R-array: (vertex, surface, subject)
    dat = R.array(
        dat, dim = [ len(hdr), len(vtx), len(sbj)],
        dimnames = [ hdr, vtx, sbj])

    ## permute the R-array: (subject, surface, vertex)
    dat = R.aperm(dat, [2,1,3])
    R.saveRDS(dat, fo)

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
    nbs = wrk['nbs']     # neighbor table for each vertex
    sz = wrk['sz']       # region size

    ## lists of vertex index neighborhoods, and output paths
    lfo, lvi = [], []
    for hm, nb, cv in izip(hms, nbs, cvs):
        lvi.append(__wm_neighbor__(nb, cv, sz))
        lfo.append(pt.join(dst, '{}{:05X}'.format(hm, cv)))
        
    ## lists of surfaces to be sampled, and subjects
    lsf, lsb = [[] for i in xrange(len(lvi))], []

    ## iterate all subjects
    print 'xt: sample ', len(lvi), 'WM areas from ', src, ':'
    for fn in gg(pt.join(src, '*')):
        if not fn.endswith('wm.npz'):
            continue
        sb = pt.basename(fn).split('.')[0]
        lsb.append(sb)
        print sb
        wm = np.load(fn)

        ## sample surfaces for subject {sb}
        ## si: surface index, hm: hemisphere, vi: vertex indices
        for si, hm, vi in izip(xrange(len(lvi)), hms, lvi):
            lsf[si].append(wm[hm][vi])

    ## write the samples to file in numpy format.
    print 'xt: write WM samples to ', dst, ':'
    sbj = np.array(lsb)
    for sf, vi, fo in izip(lsf, lvi, lfo):
        np.savez_compressed(fo + '.npz', sbj=sbj, vtx=np.array(sf))
        vi = ['{:05X}'.format(i) for i in vi]
        __sample2rds__(sf, lsb, vi, fo + '.rds')
        print fo + ": created"

    print 'xt: success'

if __name__ == "__main__":
    pass
    import sys
    if len(sys.argv) < 2:
        wrk = '../tmp/wm_samples/0B_0078/script/000_00.pk'
        with open(wrk) as pk:
            wrk = cPickle.load(pk)
    else:
        wrk = sys.argv[1]
        with open(wrk) as pk:
            wrk = cPickle.load(pk)
        __sample_wm__(wrk)
