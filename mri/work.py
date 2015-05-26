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
        action = 'created'
        if pt.isfile(fo):
            if not ovr:
                print fo, "exists"
                continue
            else:
                action = 'renewed'

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

        print fo, action

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
    
def vtx2vox(src, dst, ovr = False, sz = 1):
    """ vertex into grid """
    hlp.mk_dir(dst)
    print "vtx2grd: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        ## offset to 0, get voxel position
        ps = sf['pos']          # reference, not a copy!
        for a in 'xyz':
            ps[a] = np.rint((ps[a] - ps[a].min()) / sz)

        ## floating point to integer position
        sf = np.array(sf, dtype = VOX)

        ## sort by position
        sf = sf[sf['pos'].argsort()]
        
        ## get unique position and starting indices
        U, S = np.unique(sf['pos'], return_index = True)

        ## combin vertices fall into the same voltex
        ## container to hold the combined vertices
        C = np.empty(U.shape, sf.dtype)   

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
            
        with open(fo, 'wb') as pk:
            cPickle.dump(C, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"    

def sfr2vlm(src, dst, ovr = False, dim = DIM):
    hlp.mk_dir(dst)

    print "srf2vlm: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        ## non-zero position
        pos = [sf['pos'][a] for a in 'xyz']
        vlm = np.zeros(dim, dtype = VLM)    # 3D volumn
        for f in VLM.names:
            vlm[f][pos]=sf[f]

        ## value of the surface
        val = fsf(sf)

        sf = {'vlm':vlm, 'val':val}
   
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)
        
        if renew:
            print fo, "renewed"
        else:
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

def extract_region(src, lbl, val):
    root = pt.join(pt.dirname(src), str(lbl))
    npy = pt.join(root, 'npy')
    grd = pt.join(root, 'grd')
    srt = pt.join(root, 'srt')
    cmb = pt.join(root, 'cmb')
    vlm = pt.join(root, 'vlm')

    vfilter(src, npy, ovr = 0, flt = lambda v: v['lbl'] == lbl)
    vtx2grd(npy, grd, ovr = 0, sz = 1)
    srt_pos(grd, srt, ovr = 0)
    cmb_pos(srt, cmb, ovr = 0)
    sfr2vlm(cmb, vlm, ovr = 0, svl = val)

    
def test():
    csv2npy('dat/csv', 'dat/npy', ovr = 0)
    vfilter('dat/npy', 'dat/tmp', ovr = 0, flt = lambda v: v['lbl'] == 2003)
    vtx2vox('dat/tmp', 'dat/grd', ovr = 1)
    # extract_region('dat/npy/*', 1003, 0) 
    # extract_region('dat/npy/*', 1035, 2) 
    # extract_region('dat/npy/*', 2003, 3) 
    # extract_region('dat/npy/*', 2035, 5)

    # vlm2trn('dat/1003/vlm', 'dat/1003/trn')

if __name__ == "__main__":
    test()
        
