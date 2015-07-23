import pdb
import numpy as np
import sys
import os
import os.path as pt
from glob import glob as gg
import cPickle
from itertools import izip
import rpy2
import rpy2.robjects as robjs
import rpy2.rlike.container as rlc
R = rpy2.robjects.r

from rpy2.robjects.packages import importr
mtxr = importr('Matrix')

def __neighbor__(nb, cv, sz):
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
            ## remove extra vertices from last round
            idx = idx[:sz]                
            mrk.intersection_update(idx)
            break

    ## allocate connectivity matrix
    cmx = np.zeros((sz, sz), dtype = '<u1')

    ## the mapping from globe index to sampled region local index
    imp = dict(zip(idx, xrange(sz)))

    ## connection tuples
    cnn = [(imp[i], imp[j]) for i in idx for j in nb[i][2:] if j in mrk]

    ## list of 2-tuples to tuple of 2-lists
    cnn = zip(*cnn)
    cmx[cnn]=1

    ## return vertex indices and connection matrix
    return (idx, cmx)
    

def __save_rds__(sfs, sid, vid, cmx, fo):
    """
    save surface samples to R binary.
    sfs: samples of a surface region.
    sid: subject IDs, which indices the 1st dimension of {sfs}
    vid: vertex IDs
    cmx: connection matrix for vertices
    """
    hdr = list(sfs[0].dtype.names)

    ## flatten the data: subject * vertex * feature
    sfs = [f for s in sfs for v in s for f in v.item()]
    sfs = R.unlist(sfs)

    ## dimension names
    names = rlc.TaggedList([hdr, vid, sid], tags = ['ftr', 'vtx', 'sbj'])
    
    ## create 3D R-array: (vertex, surface, subject)
    sfs = R.array(
        sfs, dim = [len(hdr), len(vid), len(sid)],
        dimnames = names)

    ## the vertex connection matrix
    names = rlc.TaggedList([vid, vid], tags = ['v.i', 'v.j'])
    cmx = robjs.IntVector([r.item() for r in cmx.flat])
    cmx = mtxr.Matrix(
        cmx, nrow = len(vid), ncol = len(vid),
        sparse = 1,
        dimnames = names)

    ## permute the R-array: (subject, surface, vertex)
    ## sfs = R.aperm(sfs, [2,1,3])
    ret = rlc.TaggedList([sfs, cmx], tags = ['sfs', 'cmx'])
    R.saveRDS(ret, fo)

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

    ## lists of output path, vertex indices, and connetion matrices
    lfo, lvi, lcn = [], [], []
    for hm, nb, cv in izip(hms, nbs, cvs):
        vi, cn = __neighbor__(nb, cv, sz)
        lvi.append(vi)
        lcn.append(cn)
        lfo.append(pt.join(dst, '{}{:05X}'.format(hm, cv)))
        
    ## lists of surfaces to be sampled, and subjects
    lsf, lsb = [[] for i in xrange(len(lvi))], []

    ## iterate all subjects
    print 'xt: sample ', len(lvi), 'WM areas from ', src, ':'
    for fn in gg(pt.join(src, '*')):
        if not fn.endswith('npz'):
            continue
        sb = pt.basename(fn).split('.')[0]
        lsb.append(sb)
        print sb
        sys.stdout.flush()
        wm = np.load(fn)

        ## sample surfaces for subject {sb}
        ## si: surface index, hm: hemisphere, vi: vertex indices
        for si, hm, vi in izip(xrange(len(lvi)), hms, lvi):
            lsf[si].append(wm[hm][vi])

    ## write the samples to file in numpy format.
    print 'xt: write WM samples to ', dst, ':'
    sys.stdout.flush()
    sbj = np.array(lsb)
    for sf, vi, cn, fo in izip(lsf, lvi, lcn, lfo):
        np.savez_compressed(fo + '.npz', sbj=sbj, vtx=np.array(sf), cmx=cn)
        vi = ['{:05X}'.format(i) for i in vi]
        __save_rds__(sf, lsb, vi, cn, fo + '.rds')
        print fo + ": created"
        sys.stdout.flush()

    print 'xt: success'
    sys.stdout.flush()


if __name__ == "__main__":
    pass
    import sys
    if len(sys.argv) < 2:
        wrk = '/tmp/WMS_0032.ppk'
        with open(wrk) as pk:
            wrk = cPickle.load(pk)
    else:
        wrk = sys.argv[1]
        with open(wrk) as pk:
            wrk = cPickle.load(pk)
        __sample_wm__(wrk)
