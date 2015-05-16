## brain surface
import csv
import os
import cPickle

## vertex fields in the surface CSV table
## (index, caption, parser)
VFD = [
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


class Surface(list):
    """ brain surface

    An list of vertices sorted by position
    sn: sample 64bit integer serial number
    fp: opened surface table in csv format
    cache: True to try cached object, if the object was
    not cached, a new Surface will be created. Use False
    to re-new or create the cache
    """
    
    def __init__(self, sn, fp):
        """ initializer
        fp: file opened on the surface vertex csv file
        """
        self.sn = sn

        fcsv = csv.reader(fp)
        fcsv.next()             # skip header

        ## iterate through rest of the csv, create one vertex per line;
        for line in fcsv:
            vtx = tuple([t(line[i]) for i, c, t in VFD])
            self.append(vtx)  # new Vertex

    def __str__(self):
        return self[0:5].__str__()

    def __repr__(self):
        return self[0:5].__repr__()
        
def save_pk(sf, ds = 'pck'):
    """ save a surface in python pickle format to destiniation
    directory, using the surface serial number as the file name

    sf: the surface to be saved
    ds: the destination directory, the default is 'pck/'

    the file name will be:
        {ds}/{sf.sn}
    where {sf.sn} is the serial number of sample surface
    """
    if not os.path.exists(ds):
        os.mkdir(ds)

    with open(os.path.join(ds, sf.sn), 'wb') as pk:
        cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

def load_pk(sn, sc = 'pck'):
    """ load a surface in python pickle form from source directory,
    given that the file name identical to the surface serial number

    sn: the serial number of sample surface
    sc: the source directory, the default is 'pck/'

    return loaded surface
    """

    fi = os.path.join(sc, sn)
    if not os.path.isfile(fi):
        return None
    
    with open(fi, 'rb') as fi:
        sf = cPickle.load(fi)
    return sf


def from_csv(fn):
    """ wrapper function to create Surface by specifying csv file,
    serial number can be inferred from the file name.

    return the created Surface object,
    fn: the csv file name.
    """
    sn = os.path.basename(fn)
    sn = os.path.splitext(sn)[0]
    fi = open(fn)
    sf = Surface(sn, fi)
    fi.close()
    return sf

def test():
    import time
    start = time.time()

    with open('s01.csv') as f:
        s = Surface('01', f)

    with open('s0L.csv') as f:
        s = Surface('0L', f)

    with open('s0R.csv') as f:
        s = Surface('0R', f)

    # s01 = load_pk('01')
    # s0l = load_pk('0L')
    # s0r = load_pk('0R')

    # print s0l
    # print s0r
    # print s01
    end = time.time()
    print end - start
    
if __name__ == "__main__":
    test()

