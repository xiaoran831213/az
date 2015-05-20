import cPickle
import glob
import os
import csv
import numpy as np

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

NPY = np.dtype([
    ('idx', '<i4'),
    ('lbl', '<i2'),
    ('slc', 'i1'),
    ('pos', '<f4', (3,)),
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
    ('pos', 'u1', (3,)),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

def csv2npy(src, dst, ovr = False, flt = None):
    """ read raw csv into surface in 2D list """
    if not os.path.exists(dst):
        os.mkdir(dst)

    print "\ncsv2sfr: ", src, " -> ", dst
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):
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

        ## put the list into numpy array:
        sf = np.array(sf, dtype = NPY)

        with open(fo, 'wb') as f:
            cPickle.dump(sf, f, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"
    
def pos2vox(src, dst, sz = 1, ovr = False):
    if not os.path.exists(dst):
        os.mkdir(dst)

    print "\posf2vox: ", src, " -> ", dst
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        with open(fi, 'rb') as pk:
            sf = cPickle.load(pk)

        for v in sf:
            v[3] = map(lambda e: int(round(e/sz)), v[3])

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def v_add(v1, v2):
    """ add up to vertices """
    return [
        min(v1[0], v2[0]),
        min(v1[1], v2[1]),
        max(v1[2], v2[2]),
        [   v1[3][0] + v2[3][0],
            v1[3][1] + v2[3][1],
            v1[3][2] + v2[3][2]],
        v1[4] + v2[4],
        v1[5] + v2[5],
        v1[6] + v2[6],
        v1[7] + v2[7],
        v1[8] + v2[8],
        v1[9] + v2[9]]

def v_sum(v_grp):
    return reduce(v_add, v_grp)

def v_avg(v_grp):
    v = v_sum(v_grp)
    n = len(v_grp)
    return [
        v[0],
        v[1], 
        v[2],
        [   v[3][0] / n,
            v[3][1] / n,
            v[3][2] / n],
        v[4] / n,
        v[5] / n,
        v[6] / n,
        v[7] / n,
        v[8] / n,
        v[9] / n]

def prt_pk(src, fs, ts = None, fv = 0, tv = None):
    """ print pickle files"""
    if ts == None:
        ts = fs + 1
    if tv == None:
        tv = fv + 5
        
    for fi in glob.glob(os.path.join(src, "*"))[fs:ts]:
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
    
def get_pk(src, si = 0, fv = 0, tv = None):
    """ get surface from pickle
    si: surface index
    fv: from vertex
    tv: to vertex
    """
    fi = glob.glob(os.path.join(src, "*"))[si]
    with open(fi, 'rb') as f:
        sf = cPickle.load(f)[fv:tv]
    return sf
        
def csv2sfr(src, dst, ovr = False, flt = None):
    """ read raw csv into surface in 2D list """
    if not os.path.exists(dst):
        os.mkdir(dst)

    import csv
    print "\ncsv2sfr: ", src, " -> ", dst
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        fi = open(fi, 'rb')
        reader = csv.reader(fi)
        reader.next()             # skip header
        sf = []
        for line in reader:
            v = [t(line[i]) for i, c, t in CSV]
            if flt and not flt(v):
                continue
            sf.append(v)
        fi.close()

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def srf2vox(src, dst, ovr = False, sz = 2):
    if not os.path.exists(dst):
        os.mkdir(dst)

    print "\nsrf2vox: ", src, " -> ", dst
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        with open(fi, 'rb') as pk:
            sf = cPickle.load(pk)

        for v in sf:
            v[3] = map(lambda e: int(round(e/sz)), v[3])

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def vox_agr(src, dst, ovr = False):
    """ aggregate vertices fall in the same voxel """
    if not os.path.exists(dst):
        os.mkdir(dst)

    xyz_key = lambda v: v[3]
    idx_key = lambda v: [v[1], v[0]]

    print "\nvox2agr: ", src, " -> ", dst
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            renew = True

        with open(fi, 'rb') as pk: # read input
            sf = cPickle.load(pk)

        ## sort by voxel so vertices is grouped
        sf.sort(key = xyz_key)
        
        ns = []                 # new surface
        v = sf.pop()            # initlal vertex
        grp = [v]               # initial voxel group
        key = v[3]              # initial voxel key
        while sf:
            v = sf.pop()
            if key == v[3]: 
                grp.append(v)
            else:               # see a new voxel
                avg = v_avg(grp)
                ns.append(v_avg(grp))
                grp = [v]
                key = v[3]
        else:                   # deal with last voxel
            ns.append(v_avg(grp))

        ns.sort(key = idx_key)
        with open(fo, 'wb') as pk:
            cPickle.dump(ns, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def get_surface_bound(s):
    max_x, max_y, max_z = (-32767,) * 3
    min_x, min_y, min_z = (+32767,) * 3
    for v in s:
        x, y, z = v[3][0], v[3][1], v[3][2]
        max_x, max_y, max_z = (
            max(x, max_x),
            max(y, max_y),
            max(z, max_z))
        min_x, min_y, min_z = (
            min(x, min_x),
            min(y, min_y),
            min(z, min_z))
    return (
        (min_x, min_y, min_z),
        (max_x, max_y, max_z))
            
def get_dataset_bound(src, m):
    """ get coordinate boundary """
    max_ = np.array((-32767,) * 3)
    min_ = np.array((+32767,) * 3)
    max_x, max_y, max_z = (-32767,) * 3
    min_x, min_y, min_z = (+32767,) * 3
    print "get_dataset_bound:", src
    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])

        with open(fi, 'rb') as pk: # read input
            sf = cPickle.load(pk)

        sb = get_surface_bound(sf)
        min_x = min(min_x, sb[0][0])
        min_y = min(min_y, sb[0][1])
        min_z = min(min_z, sb[0][2])
        max_x = max(max_x, sb[1][0])
        max_y = max(max_y, sb[1][1])
        max_z = max(max_z, sb[1][2])

        print sn, sb
    return (
        (min_x - m, min_y - m, min_z - m),
        (max_x + m, max_y + m, max_z + m))

def voffset(src, dst, offset, ovr = False):
    """ offset voxels, report top, right, deap corner
    """
    if not os.path.exists(dst):
        os.mkdir(dst)
    print "aln2org: ", src, " -> ", dst

    dx = -offset[0][0]                       # offset x
    dy = -offset[0][1]                       # offset y
    dz = -offset[0][2]                       # offset z

    for fi in glob.glob(os.path.join(src, "*")):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)
        renew = False
        if os.path.isfile(fo):     # skip exists
            if not ovr:
                print fo, "exists"
                continue
            else:
                renew = True

        with open(fi, 'rb') as pk: # read input
            sf = cPickle.load(pk)

        print get_surface_bound(sf), " -> ",
        for v in sf:
            v[3] = v[3][0] + dx, v[3][1] + dy, v[3][2] + dz
        print get_surface_bound(sf)

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

        if renew:
            print fo, "renewed"
        else:
            print fo, "created"

def test():
    csv2npy('dat/csv', 'dat/npy', ovr = 0)
#    csv2sfr("dat/csv", "dat/srf", ovr = 1, flt = lambda v: v[1] == 1026)
#    srf2vox("dat/srf", "dat/vox", ovr = 1, sz = 2)
#    vox_agr("dat/vox", "dat/agr", ovr = 1)
#    print get_dataset_bound("dat/agr")
#    voffset("dat/agr", "dat/aln", ovr = 1, offset = 0)
    
if __name__ == "__main__":
    test()
        
