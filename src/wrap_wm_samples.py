import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip
from shutil import copy as cp

def write_wmsmp_script(
        src, dst = 0,
        n = 8, sz = 11, sd = 120,
        nodes = 4, psz = 32, hpc = 1):
    """
    randomly pick WM regions across subjects in {src}.
    """
    dst = pt.join(pt.dirname(src), 'wm_samples') if dst is 0 else dst
    dst = pt.normpath(pt.join(dst, '{:02X}_{:04X}'.format(sz, sd)))
    pfx = 'script'
    hlp.mk_dir(pt.join(dst, pfx))
    cp('sample_wm.py', pt.join(dst, pfx))
    
    sd = 0 if sd is None else sd
    n, sz = 2 ** n, 2 ** sz

    ## pick hemisphere and center vertices, and divide them to
    ## multiple tasks
    """ randomly pick hemispheres and center vertices """
    import random
    random.seed(sd)

    ## choose hemispheres and cooresponding neighborhood
    hms = [('lh', 'rh')[random.randint(0, 1)] for i in xrange(0x8000)][:n]

    ## read neighborhood table
    vnb = hlp.get_pk(pt.join(src, 'wm.vtxnbr.ppk'))
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
    nsbj = len([f for f in os.listdir(src) if f.endswith('wm.npz')])

    ## write commands
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    tsbj = 0.00005              # processing time per subject
    mem = nodes * 1.0             # for hpc
    wtm = psz * nsbj * tsbj     # for hpc
    nbat = 0
    bsz = nodes * psz             # batch size
    cmd = 'python {p}/sample_wm.py {p}/{t}.pk &>{t}.log || []\n'
    for i in xrange(0, n, psz):
        j = i % bsz             # within batch index
        if j == 0:              # new batch
            f = open(pt.join(dst, pfx, fbat.format(nbat)), 'wb')
            if hpc:
                hlp.write_hpcc_header(
                    f, mem = mem, walltime = wtm, nodes = nodes)
                f.write('\n')
                
        ## next nodes
        icpu = j / psz
        f.write('## node {:02d}\n'.format(icpu))
        f.write('(\n')

        ## save the working material specification for one nodes line
        tsk = '{:03d}_{:02d}'.format(nbat, icpu)
        whr = '{}/{}.pk'.format(pt.join(dst, pfx), tsk)
        wrk = {
            'hms' : hms[i : i + psz],     # hemispheres
            'cvs' : cvs[i : i + psz],     # center vertices
            'nbs' : nbs[i : i + psz],     # neighborhood table
            'sz' : sz,                    # region size (# of vertices)
            'src' : pt.abspath(src),      # source directory
            'dst' : pt.abspath(dst)}      # target directory
        hlp.set_pk(wrk, whr)
                
        ## write command for the nodes
        f.write(cmd.format(p=pfx, t=tsk))

        ## end of one nodes line
        f.write(')&\n\n')
        
        ## end of one batch
        if (i + psz) % bsz == 0:
            f.write('wait\n')
            nbat += 1
            f.close()

    # the left over
    if not f.closed:
        f.write('wait\n')
        nbat += 1
        f.close()

    ## write submition script
    f = open(pt.join(dst, 'qsb.sh' if hpc else 'bsh.sh'), 'wb')
    f.write('#!/bin/bash\n')
    frun = 'qsub {}\n' if hpc else 'sh {}\n'
    for i in xrange(nbat):
        bat = pt.join(pfx, fbat.format(i))
        f.write(frun.format(bat))
    f.close()
    hlp.chmod_x(f)

def test():
    write_wmsmp_script('../tmp/wm_asc2npz', '../tmp/wm_samples')

if __name__ == "__main__":
    pass

