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

## iterator of filenames
def i_fns(src = "*"):
    if pt.isdir(src):
        src = pt.join(src, "*")
    for fi in gg(src):
        if not pt.isfile(fi):
            continue
        yield fi

## iterator of pickled objects    
def i_pks(src = "*", ssn = False):
    if pt.isdir(src):
        src = pt.join(src, "*")
    for fi in gg(src):
        if not pt.isfile(fi):
            continue
        sn = pt.basename(pt.splitext(fi)[0])
        with open(fi, 'rb') as pk:
            sf = cPickle.load(pk)
        if ssn:
            yield (sn, sf)
        else:
            yield sf

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
    for sn, sf in i_pks(src, ssn = True):
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
    for sn, sf in i_pks(src, ssn = True):
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
    for sn, sf in i_pks(src, ssn = True):
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

def sfr2vlm(src, dst, ovr = False, dim = 32):
    hlp.mk_dir(dst)

    dim = dict(zip("xyz", (dim,) * 3))
    print "srf2vlm: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue

        ## non-zero position within dimension
        pos = sf['pos']
        msk = [pos[e] < dim[e] for e in 'xyz']
        msk = msk[0] * msk[1] * msk[2]
        idx = [pos[msk][e] for e in 'xyz']

        vlm = np.zeros(dim.values(), dtype = VLM)    # 3D volumn
        for f in VLM.names:
            vlm[f][idx]=sf[f][msk]

        with open(fo, 'wb') as pk:
            cPickle.dump(vlm, pk, cPickle.HIGHEST_PROTOCOL)
        
        print fo, "created"

def vlm2trn(src, dst, ovr = False):
    hlp.mk_dir(dst)

    print "vlm2trn: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        val = sf['val']
        vlm = sf['vlm']
        msk = np.int8(vlm['lbl'] > 0)

        sf = {'x':msk, 'y':val}
   
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)
        
        if renew:
            print fo, "renewed"
        else:
            print fo, "created"
            
def s_prt(src, fs = 0, ts = None, fv = 0, tv = None):
    """ print pickle binary"""
    if ts == None:
        ts = fs + 1
    if tv == None:
        tv = fv + 5
        
    for fi in glob.glob(pt.join(src, "*"))[fs:ts]:
        print fi + ":"
        with open(fi, 'rb') as f:
            sf = cPickle.load(f)[fv:tv]
        for v in sf:
            v = list(v)
            v[0] = "{:5d}".format(v[0])
            v[1] = "{:4d}".format(v[1])
            v[2] = "{:3d}".format(v[2])
            v[3] = map(lambda e: "{:7.2f}".format(e), v[3])
            v[3] = ",".join(v[3])
            v[4] = "{:5.2f}".format(v[4])
            v[5] = "{:5.2f}".format(v[5])
            v[6] = "{:6.4f}".format(v[6])
            v[7] = "{:4.2e}".format(v[7])
            v[8] = "{:5.2f}".format(v[8])
            v[9] = "{:5.2f}".format(v[9])
            v = " ".join(v)
            print v
    print
    
def s_get(src, si = 0):
    """ get surface from pickle
    si: surface index
    fv: from vertex
    tv: to vertex
    """
    fi = list(i_fns(src))[si]
    with open(fi, 'rb') as f:
        print fi + ": fetched"
        sf = cPickle.load(f)
    return sf

def test():
    from time import time
    csv2npy('dat/csv', 'dat/npy', ovr = 0)
    t1 = time()
    
    vtx2vox('dat/npy', 'dat/vox/1003', ovr = 1, flt = lambda v: v['lbl'] == 1003)
    sfr2vlm('dat/vox/1003', 'dat/vlm/1003', ovr = 1, dim = 32)

    vtx2vox('dat/npy', 'dat/vox/1035', ovr = 1, flt = lambda v: v['lbl'] == 1035)
    sfr2vlm('dat/vox/1035', 'dat/vlm/1035', ovr = 1, dim = 32)
    
    vtx2vox('dat/npy', 'dat/vox/2003', ovr = 1, flt = lambda v: v['lbl'] == 2003)
    sfr2vlm('dat/vox/2003', 'dat/vlm/2003', ovr = 1, dim = 32)
    
    vtx2vox('dat/npy', 'dat/vox/2035', ovr = 1, flt = lambda v: v['lbl'] == 2035)
    sfr2vlm('dat/vox/2035', 'dat/vlm/2035', ovr = 1, dim = 32)

    t2 = time()
    print t2 - t1

if __name__ == "__main__":
    test()
        
