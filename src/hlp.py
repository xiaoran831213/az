## comman helper functions
import cPickle
from glob import glob as gg
import os.path as pt
import os
import numpy as np
import sys
import pdb

def itr_fn(src = "", fmt = 'n', flt = None, drop = True):

    """
    filename iterator

    drop: drop list structure if only one file attribute
    is returned.
    format code:
    n: file name, N: absolute file name
    c: core name, C: absolute core name
    b: base name, B: absolute base name
    d: directory, D: absolute directory
    e: extension, E: absolute extension
    """
    src = resolve_path(src)
    if pt.isdir(src):
        src = pt.join(src, "*")

    if flt == None:
        flt = lambda w: True

    i = 0
    for fn in gg(src):
        if not flt(fn):
            continue
        rt = []
        for c in fmt:
            if c == 'i':
                r = i
            elif c == 'n':
                r = fn
            elif c == 'N':                # absolute filename
                r = pt.abspath(fn)
            elif c == 'C':                # absolute corename
                r = pt.abspath(fn).split('.')[0]
            elif c == 'c':                # ralative corename
                r = pt.basename(fn).split('.')[0]
            elif c == 'B':                # basename.extension
                r = pt.basename(pt.abspath(fn))
            elif c == 'b':                # basename
                r = pt.basename(fn)       
            elif c == 'D':                # absolute directory
                r = pt.dirname(pt.abspath(fn))
            elif c == 'd':                # relative directory
                r = pt.dirname(fn)
            elif c == 'e':                # extension(s)
                r = pt.basename(fn).split('.')[1:]
                if len(r) == 1:
                    r = r[0]
                if len(r) == 0:
                    r = None
            elif c == 'E':
                r = pt.basename(fn).split('.')[1:]
                if len(r) > 0:
                    r = r[-1]
                if len(r) == 0:
                    r = None
            else:
                continue
            rt.append(r)
        i += 1
        if drop and len(rt) == 1:
            yield rt[0]
        else:
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
    for rt in itr_fn(src, fmt, drop = False):
        with open(rt[0], 'rb') as pk:
            rt[0] = cPickle.load(pk)
        yield rt

def mk_dir(d):
    """ make deep folder """
    try:
        os.makedirs(d)
    except OSError as e:
        if not e.args[1] == 'File exists':
            raise e

def AUC(x, z):
    from sklearn.metrics import roc_auc_score
    x = x.reshape(x.shape[0], -1)
    z = z.reshape(z.shape[0], -1)
    s = np.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s.mean()

def resolve_path(path, full = False):
    path = pt.expanduser(path)
    path = pt.expandvars(path)
    if full:
        path = pt.abspath(path)
    return path

def write_hpcc_header(fo = None, mem = 4, walltime = 4, nodes = 8, ppn = 1):
    if fo == None:
        fo = sys.stdout
    fo.write('#!/bin/bash -login\n')
    fo.write('#PBS -l nodes={}:ppn={}\n'.format(nodes, ppn))
    hh = int(walltime);
    walltime -= hh
    mm = int(walltime * 60);
    fo.write('#PBS -l walltime={:02d}:{:02d}:00\n'.format(hh, mm))
    fo.write('#PBS -l mem={}M\n'.format(int(mem*1024)))
    fo.write('#PBS -j oe\n')

    fo.write('cd $PBS_O_WORKDIR')
    fo.write('\n')

def chmod_x(fi):
    ## make the script executable
    import stat as st
    fn = fi.name if isinstance(fi, file) else fi
    mode = os.stat(fn).st_mode
    os.chmod(fn, mode|st.S_IXUSR|st.S_IXGRP|st.S_IXOTH)
    
if __name__ == "__main__":
    pass
