## brain surface
import csv
import os
import cPickle
import glob

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

##
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
        return self.sn + "\n" + "\n".join(v.__str__() for v in self[0:5])

    def __repr__(self):
        return self[0:5].__repr__()

##
def infer_sn(fr):

    """ infer subject serial number from given object
    fr: object from which to infer serial number, it could be
    the path to a csv surface vertex table, or the file object
    open on such a table.
    """
    
    if isinstance(fr, str):
        fn = fr
    elif isinstance(fr, file):
        fn = fr.name
    else:
        fn = str(fr)
    return os.path.splitext(os.path.basename(fn))[0]

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
def del_cache(sn = None, cd = None):
    """ delete cached Surface object
    sn: wildcard pattern of subject serial numbers, by
    default it is '*', which means all cached surfaces.

    cd: cache directory, default is 'tmp/'.
    """
    if not sn:
        sn = '*'
    if not cd:
        cd = 'tmp'
    ds = os.path.join(cd, sn)
    for f in glob.glob(ds):
        os.remove(f)

##
def test():
    import time
    import glob

    ## try parse some csv vertex tables
    start = time.time()
    
    lsf = []
    for sc in glob.glob("dat/csv/*"):
        sf = make(sc, rt = False)
        lsf.append(sf)

    print "parse csv, time = ", time.time() - start

    ## show some surfaces
    print "\n".join(lsf)

##    
if __name__ == "__main__":
    test()
