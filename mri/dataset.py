## the whole study data set
import os
import glob
import surface

class Dataset(set):
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
        else:
            for fn in glob.glob(os.path.join(tag, "*")):
                self.add(surface.infer_sn(fn))
        self.tag = tag

        if not src:
            return

        ## parse csv surface vertex tables!
        src = os.path.join(src, "*")

        ## turn object return off(rt = False) so surface.make
        ## only build missing caches for newly seen csv tables,
        ## and only the serial number is returned.
        for fn in glob.glob(src):
            self.add(surface.make(fn, cd = tag, rt = False))

    def faces(self):
        """ return a Surface instance generator
        It fetches cached surface one by one according to subject
        serial numbers listed in this dataset.
        """
        for sn in self:
            yield surface.make(sn, cd = self.tag)

    def apply(self, dst = None, OP = None, *ls, **dc):
        """ apply opertaion to surfaces one by one

        dst: destination dataset name. the processed surface will
        be cached at the destination in python pickle files.

        op: the operation to be applied to each surface, it must
        take a Surface object as its first argument

        vb: verberose switch, if on, every return value is reported

        *ls, **dc: additional list and dictionary argumemts to
        pass to op

        returns:
        A tuple is returned, one is the destination dataset object,
        another is a list of OP's return value on all surfaces.

        If dst is None, meaning that OP is calculation without
        modification of the surface data, only a vector of OP's
        return value is returned
        If all of OP's return value is None, only the destination
        dataset is returned
        If both dst OP's return value on all the surfaces is None,
        None is returned.
        """
        ## Surface object generator
        SFG = self.faces()

        val = []
        if dst:
            if OP:
                for sf in SFG:
                    val.append(OP(sf, *ls, **dc))
                    surface.save(sf, dst)
                    print sf.sn, val[-1]
            else:
                if self.tag == dst:
                    pass
                import shutil
                for sn in self:
                    p0 = os.path.join(self.tag, sn);
                    p1 = os.path.join(dst, sn);
                    shutil.copy(p0, p1)
                    print sn
            dst = Dataset(dst)
        else:
            if OP:
                for sf in SFG:
                    val.append(OP(sf, *ls, **dc))
                    print val[-1]
            else:
                pass
            dst = None

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

c = 0
def test():
    import sys
    import time
    reload(surface)

    start = time.time()
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = "dat/csv"

    if len(sys.argv) > 2:
        tag = sys.argv[2]
    else:
        tag = "dat/ds0"
    d0 = Dataset(tag, src)
    # for sf in mris.surfaces():
    #     print sf
    print time.time() - start

    ## try the suface apply
    def f(sf):
        global c
        c += 1
    r = d0.apply('dat/ds1', f)
    print r

if __name__ == "__main__":
    test()
