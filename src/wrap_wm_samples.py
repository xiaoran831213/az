import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip

def write_wmsmp_script(src, dst = 0, n = 10, sz = 9, seed = 120):
    """
    randomly pick WM regions across subjects in {src}.
    """
    dst = pt.join(pt.dirname(src), 'wm_samples') if dst is 0 else dst
    dst = pt.normpath(pt.join(dst, '{:02X}_{:04X}'.format(sz, seed)))
    pfx = 'script'
    hlp.mk_dir(pt.join(dst, pfx))
    
    seed = 0 if seed is None else seed
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

    step = 32
    tsk = '/tmp/WMS_{i:04d}.ppk'
    cmd = 'python sample_wm.py /tmp/WMS_{i:04d}.ppk &>{i:04d}.log\n'
    for fo, i in hlp.hpcc_iter(
            xrange(0, n, step), dst, npb=4, mpn=2, tpp=0.1,
            mds=['R/3.1.0'],
            lnk=['wm_sample.py'],
            debug=True):

        ## save the working material specification for one nodes line
        wrk = {
            'hms' : hms[i : i + step],     # hemispheres
            'cvs' : cvs[i : i + step],     # center vertices
            'nbs' : nbs[i : i + step],     # neighborhood table
            'sz' : sz,                    # region size (# of vertices)
            'src' : pt.abspath(src),      # source directory
            'dst' : pt.abspath(dst)}      # target directory

        hlp.set_pk(wrk, tsk.format(i=i))
        fo.write(cmd.format(i=i))

def test():
    write_wmsmp_script('../tmp/align_vtx', '../tmp/wm_samples')
    
if __name__ == "__main__":
    pass

