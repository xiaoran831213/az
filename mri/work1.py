import cPickle
import glob
import os
import csv
import numpy as np
import os.path as pt
from scipy.sparse import dok_matrix

## csv to surface
CSV = [
    (0, "Index", int),
    (1, "Label", int),
    (2, "Sulcus", int),
    (3, "coordinates", lambda w: tuple(eval(w))),
    (4, "area", float),
    (5, "mean curvature", float),
    (6, "travel depth", float),
    (7, "geodesic depth", float),
    (8, "FreeSurfer convexity", float),
    (9, "FreeSurfer thickness", float)]

F3D = np.dtype([
    ('x', '<f4'),
    ('y', '<f4'),
    ('z', '<f4')])


I3D = np.dtype([
    ('x', '<u1'),
    ('y', '<u1'),
    ('z', '<u1')])


NPY = np.dtype([
    ('idx', '<i4'),
    ('lbl', '<i2'),
    ('slc', 'i1'),
    ('pos', F3D),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

VOX = np.dtype([
    ('idx', '<i4'),
    ('lbl', '<i2'),
    ('slc', 'i1'),
    ('pos', I3D),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

BND = np.dtype([
    ('ssn', np.str, 32),
    ('min', '<f4', (3,)),
    ('max', '<f4', (3,)),
    ('len', '<f4', (3,))])

def i_fns(fdr = "", ptn = "*"):
    for f in glob.glob(pt.join(fdr, ptn)):
        yield f

def g_cnt(fdr = "", ptn = "*"):
    return len(glob.glob(pt.join(fdr, ptn)))

def i_pks(fdr = "", ptn = "*", ssn = False):
    for fi in glob.glob(pt.join(fdr, ptn)):
        sn = pt.basename(pt.splitext(fi)[0])
        with open(fi, 'rb') as pk:
            sf = cPickle.load(pk)
        if ssn:
            yield (sn, sf)
        else:
            yield sf

def i_pos(fdr = "", ptn = "*", ssn = False):
    for fi in glob.iglob(pt.join(fdr, ptn)):
        sn = pt.basename(pt.splitext(fi)[0])
        with open(fi, 'rb') as pk:
            ps = cPickle.load(pk)['pos']
        if ssn:
            yield (sn, ps)
        else:
            yield ps

def csv2npy(src, dst, ovr = False, flt = None):
    """ read raw csv into surface in 2D list """
    if not os.path.exists(dst):
        os.mkdir(dst)

    print "\ncsv2sfr: ", src, " -> ", dst
    for fi in i_fns(src):
        sn = os.path.basename(pt.splitext(fi)[0])
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        ## open vertices csv table, skip header
        with open (fi, 'rb') as f:
            vt = csv.reader(f)
            vt.next()             

            ## load vertex table into python list
            sf = []
            for line in vt:
                v = tuple([t(line[i]) for i, c, t in CSV])
                if flt == None or flt(v):
                    sf.append(v)

        sf = np.array(sf, dtype = NPY)   # list to np.array
        with open(fo, 'wb') as f:
            cPickle.dump(sf, f, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"


def g_bnd(src):
    """ get dataset coordinate bound """
    s_m = []
    for sn, ps in i_pos(src, ssn = True):
        p_min = np.array(
            (ps['x'].min(), ps['y'].min(), ps['z'].min()),
            dtype = '<f4')
        p_max = np.array(
            (ps['x'].max(), ps['y'].max(), ps['z'].max()),
            dtype = '<f4')
        s_m.append((sn, p_min, p_max, p_max - p_min))
    s_m = np.array(s_m, dtype = BND)
    return s_m
    
def vtx2grd(src, dst, sz = 1, ovr = False):
    """ vertex into grid """
    if not pt.exists(dst):
        os.mkdir(dst)

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

        sf['slc'] = sf['slc'] > -1
        sf = np.array(sf, dtype = VOX)
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def srt_pos(src, dst, ovr = False):
    if not pt.exists(dst):
        os.mkdir(dst)
        
    print "srt_pos: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        ## sort by position and save
        p = sf['pos']
        sf = sf[np.lexsort((p[:,2], p[:,1], p[:,0]))]
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def cmb_pos(src, dst, ovr = False):
    if not pt.exists(dst):
        os.mkdir(dst)
        
    print "cmb_pos: ", src, " -> ", dst
    for sn, sf in i_pks(src, ssn = True):
        fo = pt.join(dst, sn)
        renew = False
        if pt.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"    

#b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
VDM = (64, 64, 64)
VLM = np.dtype([
    ('ssn', '<i4'),
    ('lbl', '<i2', VDM),
    ('slc', '<u1', VDM),
    ('crv', '<f4', VDM),
    ('are', '<f4', VDM),
    ('tdp', '<f4', VDM),
    ('gdp', '<f4', VDM),
    ('cnv', '<f4', VDM),
    ('tck', '<f4', VDM)])

def sfr2vlm(src, dst, dim, ovr = False):
    if not pt.exists(dst):
        os.mkdir(dst)

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
        p_min = sf['pos'].min(axis = 0)

        ## offset to 0, get voxel position
        sf['slc'] = sf['slc'] > -1
        sf = np.array(sf, dtype = VOX)
        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)
        
        if renew:
            print fo, "renewed"
        else:
            print fo, "created"
    

def prt_sf(src, fs = 0, ts = None, fv = 0, tv = None):
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
    
def get_sf(src, si = 0, fv = 0, tv = None):
    """ get surface from pickle
    si: surface index
    fv: from vertex
    tv: to vertex
    """
    fi = glob.glob(pt.join(src, "*"))[si]
    with open(fi, 'rb') as f:
        print fi + ":"
        sf = cPickle.load(f)[fv:tv]
    return sf
        

def get_ps(src, si = 0):
    return get_sf(src, si)['pos']


def test():
#    csv2npy('dat/csv', 'dat/npy1', ovr = 1, flt = lambda v: v[1] == 1011)
#    vtx2grd('dat/npy1', 'dat/grd1', ovr = 1, sz = 1)
    srt_pos('dat/gr1d', 'dat/srt1', ovr = 1)
#    print get_dataset_bound("dat/agr")
#    voffset("dat/agr", "dat/aln", ovr = 1, offset = 0)
    
if __name__ == "__main__":
    test()
        
