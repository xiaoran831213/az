import cPickle
import glob

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

VTX = [
    (0, "Index", int),
    (1, "Label", int),
    (2, "Sulcus", int),
    (3, "x", float),
    (4, "y", float),
    (5, "z", float),
    (6, "area", float),
    (7, "mean curvature", float),
    (8, "travel depth", float),
    (9, "geodesic depth", float),
    (10, "FreeSurfer convexity", float),
    (11, "FreeSurfer thickness", float)]
    
def get_pk(src, idx):
    for i, f in enumerate(glob.glob(os.path.join(src, "*"))):
        if i < idx:
            continue
        with open(f, 'rb') as pk:
            return cPickle.load(pk)
    
def csv2sfr(src, dst):
    if not os.path.exists(dst):
        os.mkdir(dst)

    for fi in glob.glob(os.path.join(src)):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)

        fi = open(fi, 'rb')
        reader = csv.reader(fp)
        reader.next()             # skip header
        sf = []
        for line in reader:
            vtx = []
            for e in (t(line[i]) for i, c, t in CSV):
                if isinstance(e, list):
                    vtx.extend(e)
                else:
                    vtx.append(e)
            data.append(vtx)
        fp.close()

        with open(fo, 'wb') as pk:
            cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

def sfr2vox(src, dst, sz = 2):
    if not os.path.exists(dst):
        os.mkdir(dst)

    IX, IY, IZ = [i for i, c, t in VTX if c in 'xyz']
    for fi in glob.glob(os.path.join(src)):
        sn = os.path.basename(os.path.splitext(fi)[0])
        fo = os.path.join(dst, sn)

        fi = open(fi, 'rb')
        sf = cPickle.load(fi)
        for v in sf:
            v[IX] = int(round(v[IX]/sz))
            v[IY] = int(round(v[IY]/sz))
            v[IZ] = int(round(v[IZ]/sz))

        
