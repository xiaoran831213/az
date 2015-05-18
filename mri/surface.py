## brain surface
import csv
import os
import cPickle
import glob

## vertex fields in the surface CSV table
## (index, caption, parser)
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

##
class Surface(list):

    """ brain surface

    An list of vertices sorted by position
    sn: sample 64bit integer serial number
    fn: surface table in csv format
    cache: True to try cached object, if the object was
    not cached, a new Surface will be created. Use False
    to re-new or create the cache
    """
    
    def __init__(self, fn):
        """ initializer
        fp: file opened on the surface vertex csv file
        """
        fn = os.path.basename(fn);
        self.sn = os.path.splitext(fn)[0]

        fp = open(fn, 'rb')
        fcsv = csv.reader(fp)
        fcsv.next()             # skip header

        ## iterate through rest of the csv, create one vertex per line;
        for line in fcsv:
            vtx = []
            for e in (t(line[i]) for i, c, t in CSV):
                if isinstance(e, list):
                    vtx.extend(e)
                else:
                    vtx.append(e)
            self.append(vtx)
        fp.close()

    def __str__(self):
        return self.sn + "\n" + "\n".join(v.__str__() for v in self[0:5])

    def __repr__(self):
        return self[0:5].__repr__()

##
def make(fr, cd = None, rt = True):
    """ helper for creation of Surface object

    sn: the serial number
    fr: from what the surface object is to be created.
        a serial nmmber means only search the cache
        a path to an csv vertex table or opened csv file means
    create surface object from the file if not already cached,
    and the serial number is inferred from the file name.
    
    cd: the Cache Directory to hold python pickle files.
    It will be searched to speed up "creation" of Surface
    sbjects. The default cache directory is 'tmp/'

    tr: True to return the Surface instance newly created or
    loaded from cached.

    False to return only the serial number.
    
    If instance return is off, the function serves as a cache
    builder only creating and caching newly encountered source
    surfaces(who does not have a cache).
    """

    if not cd:
        cd = 'tmp'

    if not os.path.exists(cd):
        os.mkdir(cd)

    sn = infer_sn(fr)
    pk = os.path.join(cd, sn)          # cache filename

    ## "create" from cache if it was cached
    if os.path.isfile(pk):
        if rt:
            with open(pk, 'rb') as pk:
                return cPickle.load(pk)
        else:
            return sn

    ## create from file, first get csv file object as fi
    if isinstance(fr, str):     # open file
        fi = open(fr)
    elif type(fr) == file:      # already opened
        fi = fr
    else:
        raise IOError(fr, 'not a path or opened file')

    ## create Surface object
    sf = Surface(sn, fi)
    if isinstance(fr, str):     # close file
        fi.close()

    ## cache the Surface object.
    with open(pk, 'wb') as pk:
        cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)

    if rt:
        return sf
    else:
        return sn

##
def save(sf, cd):
    """ make cache for Surface object

    sf: the Surface object to be saved
    cd: cache directory to save the object. it will be
    automatically created.
    """
    if not os.path.exists(cd):
        os.mkdir(cd)

    pk = os.path.join(cd, sf.sn)       # cache filename

    ## cache the Surface object
    with open(pk, 'wb') as pk:
        cPickle.dump(sf, pk, cPickle.HIGHEST_PROTOCOL)
    
##
def test():
    import time
    import glob

    ## try parse some csv vertex tables
    start = time.time()

    for i, sc in enumerate(glob.glob("dat/csv/*")):
        t1 = time.time()
        print make(sc, rt = False), time.time() - t1
    print "parse csv, time = ", time.time() - start, "\n"

    print "show some parsed surface:"
    for i, sc in enumerate(glob.glob("dat/csv/*")):
        t1 = time.time()
        sf = make(sc)
        print sf.sn, "load time: ", time.time() - t1
        print sf, "\n"
    
##    
if __name__ == "__main__":
    test()
