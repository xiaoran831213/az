## comman helper functions
import cPickle
from theano import shared as ts
from glob import glob as gg
import os.path as pt
import os
import numpy as np

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

def get_ts(src, idx = 0, tag = None, dtp = None):
    """ load data into a theano.shared object. """
    if pt.isdir(src):
        fn = gg(pt.join(src, "*"))[idx]
    else:
        fn = src
        
    with open(fn, 'rb') as f:
        print fn + ": fetched"
        obj = cPickle.load(f)

    if dtp:
        obj = np.asarray(obj, dtype = dtp)
    if not tag:
        tag = fn
        
    return ts(obj, name = tag, borrow = True)

def AUC(x, z):
    from sklearn.metrics import roc_auc_score
    x = x.reshape(x.shape[0], -1)
    z = z.reshape(z.shape[0], -1)
    s = np.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s.mean()

def num_pk(src):
    return len(gg(src))

def itr_pk(src, bsn = False, ext = False):
    """ pickle file iterator """
    if pt.isdir(src):
        src = pt.join(src, "*")

    for fn in gg(src):
        if (pt.isdir(fn)):
            print fn, ": is a directory"
            continue
            
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
