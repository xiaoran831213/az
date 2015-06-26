import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip

def __resolve__sid__(ids = None):
    ## get subjects
    itr = hlp.itr_fn(
        "$SUBJECTS_DIR", 'c',
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
            ids = set(hlp.itr_fn(ids, 'c'))

    ids = set(ids)
    sbj.intersection_update(ids)
    return sbj
    
def write_align_script(ids = None, dst = "../hpc/align_vtx", cpu = 4, psz = 8, hpc = True):
    """
    ids: list of subject to call align_value.sh, None mean use all subjects
    in $SUBJECTS_DIR
    dst: output destination of vertex alignment
    cpu: number of sub processes (cpus) in each submit.
    psz: process size - how much command to run on one cpu.
    """
    import shutil
    tsk_dir = 'script'       # where to write scripts?
    hlp.mk_dir(pt.join(dst, tsk_dir))
    shutil.copy('align_vtx.sh', pt.join(dst, tsk_dir))
    dst=hlp.resolve_path(dst)

    ## read subject ids
    ids=__resolve__sid__(ids)

    ## container to hold the combined recon-all command
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    mem = cpu * 1.00            # for hpc
    wtm = psz * 0.06            # for hpc
    nbat = 0
    bsz = cpu * psz             # batch size
    cmd = tsk_dir + '/align_vtx.sh -s {sid} -d . &>> {sid}.log || [ ]\n'
    for i, s in enumerate(ids):
        j = i % bsz        # within batch index
        if j == 0:         # new batch
            f = open(pt.join(dst, tsk_dir, fbat.format(nbat)), 'wb')
            if hpc:
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
        if not (j + 1) % psz == 0:
            f.write(')&\n\n')
        f.write('wait\n')
        nbat += 1
        f.close()

    ## write submition script
    f = open(pt.join(dst, 'qsb.sh' if hpc else 'bsh.sh'), 'wb')
    f.write('#!/bin/bash\n')
    frun = 'qsub {}\n' if hpc else 'sh {}\n'
    for i in xrange(nbat):
        bat = pt.join(tsk_dir, fbat.format(i))
        f.write(frun.format(bat))
    f.close()

    import stat as st
    mode = os.stat(f.name).st_mode
    os.chmod(f.name, mode|st.S_IXUSR|st.S_IXGRP|st.S_IXOTH)

## type definitions
__F3D = np.dtype([('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
__NPY = np.dtype([
    ("xyz", __F3D),
    ("are", "<f4"),
    ("crv", "<f4"),
    ("sul", "<f4"),
    ("tck", "<f4")])
    
__ASC = [
    (0, "xyz", lambda w:tuple(eval(w))),
    (1, "are", float),
    (2, "crv", float),
    (3, "sul", float),
    (4, "tck", float)]

def __extract_neighbor_table__(dst = None, ovr = False):
    """
    extract vertex neighborhood table from the white matter surface
    of freesurfer atlas, and save them into python pickle format.
    """
    from subprocess import call
    from shutil import copy
    if dst == None:
        dst = "."
    else:
        hlp.mk_dir(dst)
    dst = pt.normpath(dst)

    fo = pt.join(dst, 'wm.vtxnbr.ppk')
    if pt.isfile(fo) and not ovr:
        print fo, ': exists'
        return
    
    output = {}
    s = pt.expandvars("$FREESURFER_HOME/subjects/fsaverage/surf")
    for h in ('lh', 'rh'):
        sfr = '{sf}/{hm}.white'.format(hm=h, sf=s)
        asc = '/tmp/{hm}.vbr.asc'.format(hm=h)
        cmd = ['mris_convert', '-v', sfr, asc]
        print cmd
        call(cmd)

        hem = []
        with open(asc) as f:
            for line in f:
                hem.append(tuple([int(w) for w in line.split()]))
        output[h] = hem
    hlp.set_pk(output, fo)

def wm_asc2npz(src, dst = None, ovr = False):
    """
    convert white matter vertex table from ascii
    to numpy array
    """
    dst = pt.join(pt.dirname(src), 'wm_asc2npz') if dst is None else dst
    hlp.mk_dir(dst)
        
    lh = hlp.itr_fn(src, flt = lambda w: w.endswith('lh.asc'))
    rh = hlp.itr_fn(src, flt = lambda w: w.endswith('rh.asc'))
    sb = hlp.itr_fn(src, 'c', flt = lambda w: w.endswith('rh.asc'))
    for s, l, r in zip(sorted(sb), sorted(lh), sorted(rh)):
        fo = pt.join(dst, s + '.wm.npz')
        if pt.isfile(fo) and not ovr:
            print fo, ": exists"
            continue

        ## left hempersphere
        ls = []
        with open (l, 'rb') as f:
            reader = csv.reader(f, delimiter = '\t')
            reader.next()  # skip header line
            for line in reader:
                v = [t(line[i]) for i, c, t in __ASC]
                ls.append(tuple(v))
            ls = np.array(ls, dtype = __NPY)

        ## right hempersphere
        rs = []
        with open (r, 'rb') as f:
            reader = csv.reader(f, delimiter = '\t')
            reader.next()  # skip header line
            for line in reader:
                v = [t(line[i]) for i, c, t in __ASC]
                rs.append(tuple(v))
            rs = np.array(rs, dtype = __NPY)

        ## the white matter surface in numpy
        np.savez_compressed(fo, lh=ls, rh=rs)
        print fo, ": created"
        
    ## extract vertex neighborhood table
    __extract_neighbor_table__(dst, ovr)


def wrap_wm_sample(
        src, dst = None,
        n = 7, sz = 10, seed = 120,
        cpu = 4, psz = 32, hpc = True):
    """
    randomly pick WM regions across subjects in {src}.
    """
    dst = pt.join(pt.dirname(src), 'wm_sample') if dst is None else dst
    dst = pt.normpath(pt.join(dst, '{:02X}_{:04X}'.format(sz, seed)))
    tsk = 'script'
    hlp.mk_dir(pt.join(dst, tsk))
    
    seed = 0 if seed is None else seed
    n, sz = 2 ** n, 2 ** sz

    ## pick hemisphere and center vertices, and divide them to
    ## multiple tasks
    """ randomly pick hemispheres and center vertices """
    import random
    random.seed(seed)

    ## choose hemispheres and cooresponding neighborhood
    hms = [('lh', 'rh')[random.randint(0, 1)] for i in xrange(0x8000)][:n]

    ## read vertex neighborhood information
    vnb = hlp.get_pk(pt.join(src, 'wm.vtxnbr.ppk'))
    nbs = (vnb[h] for h in hms)

    ## choose center vertices, expecting n < number of vertices
    cvs = {}
    for h in ('lh', 'rh'):
        idx = range(len(vnb['lh']))
        random.shuffle(idx)
        cvs[h] = idx[:n]
    cvs = [cvs[h][i] for i, h in enumerate(hms)]

    ## write commands
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    mem = cpu * 1.00            # for hpc
    wtm = psz * 0.06            # for hpc
    nbat = 0
    bsz = cpu * psz             # batch size
    cmd = tsk + '/sample_wm.py "{s}" "{d}" {hc} || []\n'
    for i in xrange(0, n, psz):
        j = i % bsz             # within batch index
        if j == 0:              # new batch
            f = open(pt.join(dst, tsk, fbat.format(nbat)), 'wb')
            if hpc:
                hlp.write_hpcc_header(
                    f, mem = mem, walltime = wtm, nodes = cpu)
                f.write('\n')
                
        ## next cpu
        icpu = j / psz
        f.write('## node {:02d}\n'.format(icpu))
        f.write('(\n')

        ## save hm_cv list for current batch_cpu
        h_c = '{}/{:03d}_{:02d}.pk'.format(pt.join(dst, tsk), nbat, icpu)
        hlp.set_pk((hms[i:i+psz], cvs[i:i+psz]), h_c)
                
        ## write command for the cpu
        f.write(cmd.format(s=pt.abspath(src), d='.', hc=h_c))

        ## end of one cpu line
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
        bat = pt.join(tsk, fbat.format(i))
        f.write(frun.format(bat))
    f.close()

    ## make the script executable
    import stat as st
    mode = os.stat(f.name).st_mode
    os.chmod(f.name, mode|st.S_IXUSR|st.S_IXGRP|st.S_IXOTH)

def test():
    write_align_script(dst = '../tmp/align_vtx')
    wm_asc2npz('../tmp/align_vtx', dst = '../tmp/wm_asc2npz')
    wm_sample('../tmp/wm_asc2npz', '../tmp/wm_sample', seed = 120)
    
if __name__ == "__main__":
    pass
