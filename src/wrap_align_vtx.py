import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import sys
import hlp

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

__CMD = hlp.resolve_path('align_vtx.sh', full = True)
__CMD = __CMD + ' -s {sid} -d {dst} &> /dev/null || [ ]\n'

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
    
def write_align_script(ids = None, dst = "align_vtx", cpu = 4, bsz = 32):
    """

    ids: list of subject to call align_value.sh, None mean use all subjects
    in $SUBJECTS_DIR
    cpu: number of sub processes (cpus) in each submit.
    bsz: batch size of each submit
    """
    hlp.mk_dir(dst)
    dst=hlp.resolve_path(dst, full = True)

    ## read subject ids
    ids=__resolve__sid__(ids)

    ## container to hold the combined recon-all command
    i_bat = 0
    i_cmd = 0
    f_bat = '{}/{:03d}.qs'
    mem = cpu * 1.0
    wtm = bsz / cpu * 0.03
    for sid in ids:
        ## new batch
        if i_cmd % bsz == 0:
            f = open(f_bat.format(dst, i_bat), 'wb')
            hlp.write_hpcc_header(f, mem = mem, walltime = wtm, nodes = cpu)
            f.write('\n')
            i_cpu = 0
            
        ## new cpu
        if i_cmd % cpu == 0:
            f.write('## node {:02d}\n'.format(i_cpu))
            f.write('(\n')
            
        ## write new command
        f.write(__CMD.format(sid=sid, dst=dst))
        i_cmd += 1

        ## end of one cpu
        if i_cmd % cpu == 0:
            f.write(')&\n\n')
            i_cpu += 1
        
        ## end of one batch
        if i_cmd % bsz == 0:
            f.write('wait\n')
            f.close()
            i_bat += 1

    ## write submitor
    if not f.closed:
        if not i_cmd % cpu == 0:
            f.write(')&\n\n')
        f.write('wait\n')
        f.close()
    f = open(pt.join(dst, 'script', 'tsk.sh'), 'wb')
    for i in xrange(i_bat):
        f_bat = pt.abspath(f_bat.format(dst, i))
        f.write('qsub {}\n'.format(f_bat))

if __name__ == "__main__":
    pass
