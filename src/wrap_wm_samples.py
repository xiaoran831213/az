import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip

def write_wmsmp_script(src, dst = None, n = 13, sz = 9, seed = None):
    """
    randomly pick WM regions across subjects in {src}.
    """
    seed = 0 if seed is None else seed
    src = pt.expandvars(src)
    dst = pt.join(pt.dirname(src), 'wm_samples') if dst is None else dst
    dst = pt.expandvars(dst)
    dst = pt.normpath(pt.join(dst, '{:02X}_{:04X}'.format(sz, seed)))
    
    n, sz = 2 ** n, 2 ** sz

    ## randomly pick hemispheres and center vertices
    import random
    random.seed(seed)

    ## sample hemispheres
    hms = [('lh', 'rh')[random.randint(0, 1)] for i in xrange(0x8000)][:n]

    ## samplecooresponding vertex neighborhoods
    vnb = hlp.get_pk(pt.join(src, 'vtx_nbr.ppk'))
    nbs = [vnb[h] for h in hms]

    ## choose center vertices, expecting n < number of vertices
    cvs = {}
    for h in ('lh', 'rh'):
        idx = range(len(vnb[h]))
        random.shuffle(idx)
        cvs[h] = idx[:n]
    cvs = [cvs[h][i] for i, h in enumerate(hms)]

    ## gather lists output file names for each center vertex, 
    ## also check and skip existing samples
    ext = set()

    for i, hm, cv in izip(xrange(n), hms, cvs):
        ## make output file name, check overwiting
        fo = pt.join(dst, '{}{:05X}.npz'.format(hm, cv))
        if pt.isfile(fo):
            print fo, ": exists"
            ext.add(i)

    hms[:] = [hms[i] for i in xrange(n) if not i in ext]
    cvs[:] = [cvs[i] for i in xrange(n) if not i in ext]
    nbs[:] = [nbs[i] for i in xrange(n) if not i in ext]
    n -= len(ext)

    ## count number of subjects
    nsbj = len([f for f in os.listdir(src) if f.endswith('.npz')])

    ## create task directory

    ## surface region to be sampled per task run.
    step = 32
    tsk = 'tsk/WMS_{i:04d}.ppk'
    cmd = 'python wm_sample.py tsk/WMS_{i:04d}.ppk &>{i:04d}.log\n'
    hlp.mk_dir(pt.dirname(pt.join(dst, tsk)))
    
    for fo, i in hlp.hpcc_iter(
            xrange(0, n, step), dst, npb=1, mpn=8, tpp=4.0,
            mds=['R/3.1.0'],
            lnk=['wm_sample.py'],
            debug=False):

        ## save the working material specification for one nodes line
        wrk = {
            'hms' : hms[i : i + step],     # hemispheres
            'cvs' : cvs[i : i + step],     # center vertices
            'nbs' : nbs[i : i + step],     # neighborhood table
            'sz' : sz,                    # region size (# of vertices)
            'src' : pt.abspath(src),      # source directory
            'dst' : pt.abspath(dst)}      # target directory
        hlp.set_pk(wrk, pt.join(dst, tsk).format(i=i))

        fo.write(cmd.format(i=i))

def test():
    write_wmsmp_script('$AZ_AVTX', '$AZ_WMSP')
    
if __name__ == "__main__":
    pass
