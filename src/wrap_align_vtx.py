import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
import stat
from glob import glob as gg
import sys
import hlp
import shutil

__F3D = np.dtype([('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
__NPY = np.dtype([
    ("xyz", __F3D),
    ("are", "<f4"),
    ("crv", "<f4"),
    ("sul", "<f4"),
    ("tck", "<f4")])
    
__TXT = [
    ("xyz", lambda w:tuple(list)),
    ("are", float),
    ("crv", float),
    ("sul", float),
    ("thk", float)]

def __resolve__sid__(ids = None):
    ## get subjects
    itr = hlp.itr_fn(
        "$SUBJECTS_DIR", 'b',
        flt = lambda f: not pt.islink(f) and pt.isdir(f))
    sbj = set(itr)
    if ids == None:
        return sbj

    ## get requested list of id
    if isinstance(ids, str):
        if pt.isfile(ids):
            with open(ids, 'rb') as f:
                ids = set(f.readlines())
        if pt.isdir(ids):
            ids = set(hlp.itr_fn(ids, 'b'))

    ids = set(ids)
    sbj.intersection_update(ids)
    return sbj
    
def write_align_script(ids = None, dst = "align_vtx", cpu = 4, psz = 8):
    """
    ids: list of subject to call align_value.sh, None mean use all subjects
    in $SUBJECTS_DIR
    dst: output destination of vertex alignment
    cpu: number of sub processes (cpus) in each submit.
    psz: process size - how much command to run on one cpu.
    """

    pbs = 'script'                        # PBS script directory
    hlp.mk_dir(pt.join(dst, pbs))
    shutil.copy('align_vtx.sh', pt.join(dst, pbs))
    dst=hlp.resolve_path(dst)

    ## read subject ids
    ids=__resolve__sid__(ids)

    ## container to hold the combined recon-all command
    fbat = '{:03d}.qs'
    mem = cpu * 1.00
    wtm = psz * 0.03
    nbat = 0
    bsz = cpu * psz
    cmd = pbs + '/align_vtx.sh -s {sid} -d . &>> {sid}.log || [ ]\n'
    for i, s in enumerate(ids):
        j = i % bsz        # within batch index
        if j == 0:         # new batch
            f = open(pt.join(dst, pbs, fbat.format(nbat)), 'wb')
            hlp.write_hpcc_header(
                f, mem = mem, walltime = wtm, nodes = cpu)
            f.write('\n')
            
        k = j % psz        # within cpu index
        if k == 0:         # new cpu
            f.write('## node {:02d}\n'.format(j / psz))
            f.write('(\n')
            
        ## write new command
        f.write(cmd.format(sid=s))

        ## end of one cpu line
        if (j + 1) % psz == 0:
            f.write(')&\n\n')
        
        ## end of one batch
        if (i + 1) % bsz == 0:
            f.write('wait\n')
            nbat += 1
            f.close()

    # the left over
    if not f.closed:
        f.write(')&\n\n')
        f.write('wait\n')
        nbat += 1
        f.close()

    ## write submition script
    f = open(pt.join(dst, 'tsk.sh'), 'wb')
    f.write('#!/bin/bash\n')
    for i in xrange(nbat):
        bat = pt.join(pbs, fbat.format(i))
        f.write('qsub {}\n'.format(bat))
    f.close()
    mode = os.stat(f.name).st_mode
    os.chmod(f.name, mode|stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH)

if __name__ == "__main__":
    pass

















