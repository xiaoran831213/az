## comman helper functions
import cPickle
from glob import glob as gg
import os.path as pt
import os

def get_pk(fdr, idx = 0):
    """ load one pickle from from a folder """
    fi = gg(pt.join(fdr, "*"))[idx]
    with open(fi, 'rb') as f:
        print fi + ": fetched"
        sf = cPickle.load(f)
    return sf

def itr_pk(fdr = "", ptn = "*", bfn = False):
    """ pickle file iterator """
    for fn in gg(pt.join(fdr, ptn)):
        with open(fn, 'rb') as pk:
            obj = cPickle.load(pk)
        if bfn:
            yield (
                basename(pt.splitext(fi)[0]),
                obj)
        else:
            yield obj

def mk_dir(d):
    """ make deep folder """
    l = []
    while d:
        l.append(d)
        d = pt.dirname(d)

    for d in reversed(l):
        if pt.isdir(d):
            continue
        os.mkdir(d)

if __name__ == "__main__":
    del os
    del gg
    del pt
    del mk_dir
    del itr_pk
    del get_pk
