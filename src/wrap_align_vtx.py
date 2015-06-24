import csv
import pdb
import numpy as np
import os
import os.path as pt
import stat
from glob import glob as gg
import hlp

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
            ids = set(hlp.itr_fn(ids, 'c'))

    ids = set(ids)
    sbj.intersection_update(ids)
    return sbj
    
def write_align_script(ids = None, dst = "../hpc/align_vtx", cpu = 4, psz = 8):
    """
    ids: list of subject to call align_value.sh, None mean use all subjects
    in $SUBJECTS_DIR
    dst: output destination of vertex alignment
    cpu: number of sub processes (cpus) in each submit.
    psz: process size - how much command to run on one cpu.
    """
    import shutil
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

## name of white matter surface vertex neighbor table.
__VNB = 'wm.vtxnbr.ppk' 

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

    fo = pt.join(dst, __VNB)
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

def vtx_asc2npz(src, dst = None, ovr = False):
    """ convert vertex table from ascii to numpy array
    """
    if dst == None:
        dst = src
    else:
        hlp.mk_dir(dst)
        
    lh = hlp.itr_fn(src, flt = lambda w: w.endswith('lh.asc'))
    rh = hlp.itr_fn(src, flt = lambda w: w.endswith('rh.asc'))
    sb = hlp.itr_fn(src, 'c', flt = lambda w: w.endswith('log'))
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

def vtx_sample(src, dst = None, n = 10, sz = 512, seed = None):
    """ randomly sample regions from wm across
    all subjects in the source directory """
    import random
    random.seed(seed)

    if dst == None:
        dst = src
    if seed == None:
        seed = 0
    sfx = 'n{:04X}z{:04X}s{:04X}'.format(n, sz, seed)
    dst = pt.join(dst, sfx)

    ## get vertex neighbor table first
    vnb = pt.join(src, __VNB)
    vnb = hlp.get_pk(vnb)

    hem = ['lh', 'rh']
    nvx = {}
    ivx = {}
    for h in hem:

        dict(zip(hem, [len(vnb[k]) for k in hem]))
        ivx[h] = random.sample(xrange(u), n) for u in nvx]        
    for f, s in hlp.itr_fn(
            src, fmt = 'nc',
            flt = lambda w: w.endswith('wm.npz')):

        ## poll center vertex
        h = hem[random.randint(0, 1)]
        wm = np.load(fi)
    return ivx
    pass
        
if __name__ == "__main__":
    pass
