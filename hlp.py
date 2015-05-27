## comman helper functions
import cPickle
from glob import glob as gg
import os.path as pt
import os

def get_pk(src, idx = 0):
    """ load one pickle from a source """
    fn = gg(src)[idx]
    with open(fn, 'rb') as f:
        print fn + ": fetched"
        sf = cPickle.load(fn)
    return sf

def num_pk(src):
    return len(gg(src))

def itr_pk(src, bsn = False, ext = False):
    """ pickle file iterator """
    if pt.isdir(src):
        src = pt.join(src, "*")

    for fn in gg(src):
        with open(fn, 'rb') as pk:
            obj = cPickle.load(pk)

        fn, ex = pt.splitext(fn)
        bn = pt.basename(fn)
        rt = [obj]
        if bsn:
            rt.append(bn)
        if ext:
            rt.append(ex)

        if len(rt) > 1:
            yield tuple(rt)
        else:
            yield rt[0]

def get_pk(src, idx = 0):
    """ get data from pickle """
    if pt.isdir(src):
        src = pt.join(src, "*")
        
    fi = gg(src)[idx]
    with open(fi, 'rb') as f:
        print fi + ": fetched"
        obj = cPickle.load(f)
    return obj

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
    pass
