## comman helper functions
import cPickle
from theano import shared as ts
from glob import glob as gg
import os.path as pt
import os
import numpy as np

def itr_fn(src = "", fmt = 'n'):
    """ filename iterator """
    if pt.isdir(src):
        src = pt.join(src, "*")

    for fn in gg(src):
        if (pt.isdir(fn)):
            print fn, ": is a directory"
            continue
        rt = []
        for c in fmt:
            if c == 'n':
                r = fn
            elif c == 'N':                # absolute filename
                r = pt.abspath(fn)
            elif c == 'C':                # absolute corename
                r = pt.splitext(pt.abspath(fn))[0]
            elif c == 'c':                # ralative corename
                r = pt.splitext(fn)[0]
            elif c == 'B':                # basename.extension
                r = pt.basename(fn)       
            elif c == 'b':                # basename
                r = pt.splitext(pt.basename(fn))[0]
            elif c == 'D':                # absolute directory
                r = pt.dirname(pt.abspath(fn))
            elif c == 'd':                # relative directory
                r = pt.dirname(fn)
            elif c == 'e':                # extension
                r = pt.splitext(fn)[1]
            else:
                continue
            rt.append(r)
        yield rt

def get_pk(src, idx = 0):
    """ get data from pickle """
    if pt.isdir(src):
        fn = gg(pt.join(src, "*"))[idx]
    else:
        fn = src
    
    with open(fn, 'rb') as fp:
        obj = cPickle.load(fp)

    print fn + ": fetched"
    return obj

def set_pk(obj, dst):
    if pt.isdir(dst):
        mk_dir(dst)
        fn = pt.join(dst, obj.__name__)
    else:
        mk_dir(pt.dirname(dst))
        fn = dst

    with open(dst, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
    print fn, ": dumpped"

def num_pk(src):
    return len(gg(src))

def itr_pk(src, fmt = ''):
    """ pickle file iterator """
    fmt = 'n' + fmt
    for rt in itr_fn(src, fmt):
        with open(rt[0], 'rb') as pk:
            rt[0] = cPickle.load(pk)
        yield rt

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

def AUC(x, z):
    from sklearn.metrics import roc_auc_score
    x = x.reshape(x.shape[0], -1)
    z = z.reshape(z.shape[0], -1)
    s = np.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s.mean()

if __name__ == "__main__":
    pass
