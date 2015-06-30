import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip
from shutil import copy as cp

def write_train_sda(src, zsd, dst = 0, nodes = 4, ppn = 4, psz = 4, hpc = 1):
    """
    train SDA with WM samples in {src}
    zsd: region size and random seed used to sample WM, decide the folder name
    of the samples, and the default destination folder for SDA
    """
    dst = pt.join(pt.dirname(src), 'trained_sda') if dst is 0 else dst
    xpt = pt.join(pt.dirname(src), 'encoded_wms')
        
    seed = int(zsd.split('_')[1], 16)
    src = pt.join(src, zsd)               # source
    dst = pt.join(dst, zsd)               # target
    xpt = pt.join(xpt, zsd)               # export

    pfx = 'script'
    hlp.mk_dir(pt.join(dst, pfx))
    cp('rdm/sda.py', pt.join(dst, pfx))
    cp('rdm/da.py', pt.join(dst, pfx))
    cp('rdm/t_hlp.py', pt.join(dst, pfx))
    hlp.mk_dir(xpt)

    ## gather WM surface samples, also check existing output
    sfs = []
    for sf in hlp.itr_fn(src, 'c', lambda w: w.endswith('npz')):
        fo = pt.join(dst, sf + '.rdm')
        if pt.isfile(fo):
            print fo, ': exists'
        else:
            sfs.append(sf)

    ## write commands
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    mem = nodes * 1.0             # for hpc
    wtm = psz * 0.2               # for hpc
    nbat = 0
    bsz = nodes * psz             # batch size
    if hpc:
        cmd = 'python {p}/sda.py {p}/{t}.pk 2>&1 1>{t}.log\n'
    else:
        cmd = 'python {p}/sda.py {p}/{t}.pk 2>&1 | tee {t}.log\n'

    for i, sf in enumerate(sfs):
        j = i % bsz             # within batch index
        if j == 0:              # new batch
            f = open(pt.join(dst, pfx, fbat.format(nbat)), 'wb')
            if hpc:
                hlp.write_hpcc_header(
                    f, mem = mem, walltime = wtm, nodes = nodes, ppn = ppn)
                f.write('#PBS -l feature=intel14\n')
                f.write('#PBS -j oe\n')
                f.write('module load NumPy\n')
                f.write('\n')
                f.write('export MKL_NUM_THREADS={}\n\n'.format(ppn))
                
        ## new node 
        if j % psz is 0 and nodes > 1:
            f.write('## node {:02d}\n'.format(j/psz))
            f.write('(\n')

        ## save the working material specification for one processor
        tsk = '{}/{}.pk'.format(pfx, sf)
        whr = pt.join(dst, tsk)
        ## src: wm sample directory, dst: distination directory
        ## wms: wm sample id
        wrk = {
            'src' : pt.abspath(src),
            'dst' : pt.abspath(dst),
            'xpt' : pt.abspath(xpt),
            'wms' : sf,
            'seed' : seed} 
        hlp.set_pk(wrk, whr)
                
        ## write command for one processor
        f.write(cmd.format(p=pfx, t=sf))

        ## end of the node 
        if (j + 1) % psz is 0 and nodes > 0:
            f.write(')&\n\n')
        
        ## end of one batch
        if (i + 1) % bsz == 0:
            if nodes > 1:
                f.write('wait\n')
            nbat += 1
            f.close()

    # the left over
    if not f.closed:
        if nodes > 1:
            if (j + 1) % psz is not 0:
                f.write(')&\n\n')
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
    write_train_sda('../hpc/wm_samples', '09_0078')
    pass

if __name__ == '__main__':
    pass
