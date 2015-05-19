import cPickle
import glob
import os

## csv to surface
CSV = [
    (0, "Index", int),
    (1, "Label", int),
    (2, "Sulcus", int),
    (3, "coordinates", eval),
    (4, "area", float),
    (5, "mean curvature", float),
    (6, "travel depth", float),
    (7, "geodesic depth", float),
    (8, "FreeSurfer convexity", float),
    (9, "FreeSurfer thickness", float)]

FLT = lambda v: v[1] > 0 or v[2] > 0

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

def get_pk(src, f1, f2, l1=0, l2=5):
    """ print pickle files"""
    for i, f in enumerate(glob.glob(os.path.join(src, "*"))):
        if i < f1:
            continue
        if i >= f2:
            break

        print f + ":"
        pk = open(f, 'rb')
        sf = cPickle.load(pk)[l1:l2]
        for v in sf:
            v[0] = "{:4d}".format(v[0])
            v[1] = "{:4d}".format(v[1])
            v[2] = "{:3d}".format(v[2])
            v[3] = map(lambda e: "{:7.3f}".format(e), v[3])
            v[3] = ",".join(v[3])
            v[4] = "{:5.2f}".format(v[4])
            v[5] = "{:5.2f}".format(v[5])
            v[6] = "{:6.4f}".format(v[6])
            v[7] = "{:.2e}".format(v[7])
            v[8] = "{:5.2f}".format(v[8])
            v[9] = "{:5.2f}".format(v[9])
            v = " ".join(v)
            print v
        pk.close()
    print
    
def csv2sfr(src, dst, ovr = False, flt = FLT):
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
            if not flt(v):
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
            pos = v[3]
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
            
    
def get_dataset_bound(d):
    """ get coordinate boundary """
    max_x, max_y, max_z = (-32767,) * 3
    min_x, min_y, min_z = (+32767,) * 3
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

    return (
        (min_x, min_y, min_z),
        (max_x, max_y, max_z))



def test():
#    csv2sfr("dat/csv", "dat/srf", ovr = 1, flt = lambda v: v[1] == 1026)
#    srf2vox("dat/srf", "dat/vox", ovr = 1, sz = 4)
    vox_agr("dat/vox", "dat/agr", ovr = 0)
    
if __name__ == "__main__":
    test()
        
