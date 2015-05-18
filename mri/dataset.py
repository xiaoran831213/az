## the whole study data set
import os
import glob
import surface

class Dataset(list):
    """ all MRI scaned brain surfaces in one set, backed by python
    pickle file cache.

    The dataset act like a python set of Surface object instances.
    Cached surface will be loaded into memory from python pickles
    upon request.
    
    Uncached Surface will be created on the fly and cached in tag
    directory with support from module 'surface'

    Only the serial numbers is actually stored in the dataset.
    """

    def __init__(self, tag, src = None):
        """ initializer
        Load surface vertex csv from source directory.

        tag: name tag of the dataset, cached file will be stored
        under the directory of this name.

        src: the source directoy containing surface csv files to
        be parsed.
        """

        if not os.path.isdir(tag):
            os.mkdir(tag)
        self.tag = tag

        ## request parsing csv surface tables
        if src:
            src = os.path.join(src, "*")

            ## turn instantiation  off(rt = False) so surface.make
            ## only build missing caches for new surface csv table
            for fn in glob.glob(src):
                surface.make(fn, cd = tag, rt = False)
                print 'from ', fn

        self.__update__()

    def __update__(self):
        self[:] = []
        path = glob.glob(os.path.join(self.tag, "*"))
        self.extend(surface.infer_sn(f) for f in path)
        
    def faces(self):
        """ return a Surface instance generator

        It fetches cached surface one by one according to subject
        serial numbers listed in this dataset.
        """
        for sn in self:
            yield surface.make(sn, cd = self.tag)

    ## fetch surface(s)
    def fetch(self, fr = 0, to = None):
        """ fetch instiated surfaces

        fr: from where, either a index or a serial number
        to: where to sotp fetching
        
        if a string serial number is passed to 'fr', the surface
        with that sn will be the start of fatching.
        
        if a integer is passed to 'fr', the surfaces with indeices
        >= {fr} and < {to} is returned; if {to} is unspecified,
        only one surface at the {fr} is returned
        """

        if isinstance(fr, str):
            fr = self.index(fr)

        if not to:
            to = fr + 1

        ret = []
        for sn in self[fr:to]:
            ret.append(surface.make(sn, cd = self.tag))

        if len(ret) > 1:
            return ret
        elif len(ret) > 0:
            return ret[0]
        else:
            return None


## make dataset
def make(tag, csv = None):
    """ helper to create or load surface dataset

    tag: name tag of the dataset, cached file will be stored
    under the directory of this name.

    csv: the source directoy containing surface csv files to
    be parsed.
    """
    return Dataset(tag, csv)


## surface-wise applier
def apply(OP, src, dst = None, *ls, **dc):
    """ apply opertaion to surfaces

    OP: the operation to be applied to source surfaces, it must
    take a surface as its first argument

    src: source dataset
    dst: target dataset where the processed surface will be saved.

    *ls, **dc: additional argumemts to pass to OP

    returns:
    A tuple of 1) the destination dataset object and 2) a list of
    OP return value on all surfaces.

    If {dst} is None, only the list of OP returns is returned;
    If all of OP returns is None, only the destination dataset
    is returned;
    If both of the above are None, None is returned.
    """
    val = []
    if isinstance(dst, Dataset):
        for sf in src.faces():
            val.append(OP(sf, *ls, **dc))
            surface.save(sf, dst.tag)
            print sf.sn, val[-1]
        dst.__update__()
    else:
        for sf in src.faces():
            val.append(OP(sf, *ls, **dc))
            print sf.sn, val[-1]

    if dst:
        if any(val):
            ret = (dst, val)
        else:
            ret = dst
    else:
        if any(val):
            ret = val
        else:
            ret = None
    return ret

def save(fr, to):
    """ helper to save surface dataset
    fr: the dataset to be saved
    
    to: where to save the dataset?
    """
    import shutil
    for sn in fr:
        p0 = os.path.join(fr.tag, sn);
        p1 = os.path.join(to.tag, sn);
        shutil.copy(p0, p1)
    return Dataset(to.tag)      # update surface list
    
def test():
    import sys
    import time
    reload(surface)

    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = "dat/csv"

    if len(sys.argv) > 2:
        tag = sys.argv[2]
    else:
        tag = "dat/ds0"

    start = time.time()
    d0 = make(tag, src)
    print "load csv to ds0, time = ", time.time() - start

    ## try the suface apply
    start = time.time()
    d1 = make("dat/ds1")
    f = lambda s: len(s)
    print apply(lambda s: len(s), src = d0, dst = d1)
    print "proc ds0 to ds1, time = ", time.time() - start, "\n"

    ## try load the destination dataset
    print 'load desination dataset'

    print d1, "\n"

    ## try get individual surface from dataset
    print 'fetch the first surface'
    print d1.fetch(0), "\n"

if __name__ == "__main__":
    test()
