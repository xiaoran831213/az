import cPickle
from glob import glob as gg
import os
import numpy as np
import os.path as pt
import sys
if not sys.path.count(pt.abspath('..')):
    sys.path.insert(0, pt.abspath('..'))
import hlp
from defs import *
import pdb

## iterator of filenames

def csv2npy(src, dst, ovr = False):
    """ read raw csv into surface in 2D list """
    import csv
    hlp.mk_dir(dst)

    print "\ncsv2sfr: ", src, " -> ", dst
    for fi in i_fns(src):
        sn = os.path.basename(pt.splitext(fi)[0])
        fo = pt.join(dst, sn)
        if pt.isfile(fo) and not ovr:
                print fo, "exists"
                continue

        ## open vertices csv table, skip header
        with open (fi, 'rb') as f:
            vt = csv.reader(f)
            vt.next()  # skip header line

            ## load vertex table into python list
            sf = []
            for line in vt:
                v = tuple([t(line[i]) for i, c, t in CSV])
                sf.append(v)

        sf = np.array(sf, dtype = NPY)   # list to np.array
        with open(fo, 'wb') as f:
            cPickle.dump(sf, f, cPickle.HIGHEST_PROTOCOL)

        print fo, "created"

def vfilter(src, dst, flt, ovr = False):
    hlp.mk_dir(dst)
    print "vfilter: ", src, " -> ", dst
    for sf, sn in hlp.itr_pk(src, fmt = 'b'):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        msk = flt(sf)
        if isinstance(msk, np.ndarray):
            sf = sf[msk]
        else:
            sf = np.empty(0, dtype = NPY)
            
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"
    
def g_bnd(src):
    """ get dataset coordinate bound """
    s_m = []
    for sf, sn in hlp.itr_pk(src, fmt = 'b'):
        ps = sf['pos']
        p_min = np.array(
            (ps['x'].min(), ps['y'].min(), ps['z'].min()),
            dtype = '<f4')
        p_max = np.array(
            (ps['x'].max(), ps['y'].max(), ps['z'].max()),
            dtype = '<f4')
        s_m.append((sn, p_min, p_max, p_max - p_min))
    s_m = np.array(s_m, dtype = BND)
    return s_m
    
def vtx2vox(src, dst, ovr = False, sz = 1, flt = None):
    """ vertex into grid """
    hlp.mk_dir(dst)

    print "vtx2grd: ", src, " -> ", dst
    for sf, sn in hlp.itr_pk(src, fmt = 'b'):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo) and not ovr:
            print fo, "exists"
            continue

        ## apply filter:
        if flt:
            sf = sf[flt(sf)]
            
        ## offset to 0, get voxel position
        ps = sf['pos']          # reference, not a copy!
        for a in 'xyz':
            ps[a] = np.rint((ps[a] - ps[a].min()) / sz)

        ## floating point to integer position
        sf = np.array(sf, dtype = VOX)

        ## combin vertices fall into the same voltex
        ## 1) sort by position
        sf = sf[sf['pos'].argsort()]
        
        ## get unique position and starting indices
        U, S = np.unique(sf['pos'], return_index = True)

        ## container to hold the combined vertices
        C = np.empty(U.shape, sf.dtype)   

        g = np.split(sf, S[1:])
        for k, g in enumerate(np.split(sf, S[1:])):
            c = C[k]
            c['idx'] = g['idx'].min();
            c['lbl'] = g['lbl'].max();
            
            ## these features are averaged
            for f in VLM.names[2:]:
                c[f] = g[f].mean()

            ## any vertex in Sulcus mean the whole group is in
            c['slc'] = any(g['slc'] > -1)

            ## and the position of the group is shared
            c['pos'] = U[k]
        else:
            sf = C
            
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        print fo, "created"    

def sfr2vlm(src, dst, ovr = False, dim = 32, lth = 0.05):
    hlp.mk_dir(dst)

    dim = dict(zip("xyz", (dim,) * 3))
    print "srf2vlm: ", src, " -> ", dst
    for sf, sn in hlp.itr_pk(src, fmt = 'b'):
        fo = pt.join(dst, sn)
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue

        ## non-zero position within dimension
        pos = sf['pos']
        msk = [pos[e] < dim[e] for e in 'xyz']
        msk = msk[0] * msk[1] * msk[2]
        n0 = msk.shape[0]
        n1 = np.count_nonzero(msk)
        lpt = (n0 - n1) / np.float32(n0)
        print "{:4.3f}".format(lpt),
        idx = [pos[msk][e] for e in 'xyz']

        vlm = np.zeros(dim.values(), dtype = VLM)    # 3D volumn
        for f in VLM.names:
            vlm[f][idx]=sf[f][msk]

        with open(fo, 'wb') as pk:
            cPickle.dump(vlm, pk, cPickle.HIGHEST_PROTOCOL)
        
        print fo, "created"

def vlmCrop(src, dim, vbs = 1):
    """ calculate voxel loss due to dimension cropping """
    dim = dict(zip("xyz", (dim,) * 3))
    if vbs > 0:
        print "vlmLoss: ", src,
    crp = np.float32(0.0)                  # sum of loss
    ssz = np.float32(0.0)
    if vbs > 1:
        print "\n",
    else:
        print " ",
    for sf, fn in hlp.itr_pk(src, fmt = 'n'):
        ## non-zero position within dimension
        pos = sf['pos']
        msk = [pos[e] < dim[e] for e in 'xyz']
        msk = msk[0] * msk[1] * msk[2]
        n0 = msk.shape[0]
        n1 = np.count_nonzero(msk)
        cr = (n0 - n1) / np.float32(n0)   # loss
        if vbs > 1:
            print "{0:s}:{1:4.3f}".format(fn, cr)
        crp += cr
        ssz += 1

    avg = crp / ssz
    if vbs > 0:
        print avg
    else:
        print "\n",
    return avg

def vlm2vmk(src, dst, ovr = False):
    """ voxels to surface masks """
    hlp.mk_dir(dst)
    print "vlm2smk: ", src, " -> ", dst

    for vlm, ssn in hlp.itr_pk(src, fmt = 'b'):
        fo = pt.join(dst, ssn)
        if pt.isfile(fo) and not ovr:
            print fo, "exists"
            continue

        msk = np.uint8(vlm['lbl'] > 0)
        with open(fo, 'wb') as pk:
            cPickle.dump(msk, pk, cPickle.HIGHEST_PROTOCOL)
        print fo, "created"

def pack(src, dst, ovr):
    hlp.mk_dir(pt.dirname(dst))
    if pt.isfile(dst) and not ovr:
        print dst, 'exists'
        return

    print "pack: ", src, " -> ", dst
    pck = [dat  for dat in hlp.itr_pk(src)]
    pck = np.array(pck)
    
    with open(dst, 'wb') as pk:
        cPickle.dump(pck, pk, cPickle.HIGHEST_PROTOCOL)

    print dst, "created"
            
def test():
    from time import time
    # csv2npy('dat/csv', 'dat/npy', ovr = 0)
    t1 = time()
    
    vtx2vox('dat/npy', 'dat/vox/1003', ovr = 1, flt = lambda v: v['lbl'] == 1003)
    sfr2vlm('dat/vox/1003', 'dat/vlm/1003', ovr = 1, dim = 48)
    # vlm2vmk('dat/vlm/1003', 'dat/vmk/1003', ovr = 1)
    # pack('dat/vmk/1003', 'dat/pck/1003', ovr = 1)

    vtx2vox('dat/npy', 'dat/vox/1035', ovr = 1, flt = lambda v: v['lbl'] == 1035)
    sfr2vlm('dat/vox/1035', 'dat/vlm/1035', ovr = 1, dim = 48)
    # vlm2vmk('dat/vlm/1035', 'dat/vmk/1035', ovr = 1)
    # pack('dat/vmk/1035', 'dat/pck/1035', ovr = 1)
    
    vtx2vox('dat/npy', 'dat/vox/2003', ovr = 1, flt = lambda v: v['lbl'] == 2003)
    sfr2vlm('dat/vox/2003', 'dat/vlm/2003', ovr = 1, dim = 48)
    # vlm2vmk('dat/vlm/2003', 'dat/vmk/2003', ovr = 1)
    # pack('dat/vmk/2003', 'dat/pck/2003', ovr = 1)
    
    vtx2vox('dat/npy', 'dat/vox/2035', ovr = 1, flt = lambda v: v['lbl'] == 2035)
    sfr2vlm('dat/vox/2035', 'dat/vlm/2035', ovr = 1, dim = 48)
    # vlm2vmk('dat/vlm/2035', 'dat/vmk/2035', ovr = 1)
    # pack('dat/vmk/2035', 'dat/pck/2035', ovr = 1)

    t2 = time()
    print t2 - t1

if __name__ == "__main__":
    pass
    #test()
        
