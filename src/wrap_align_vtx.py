import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp

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
    wtm = psz * 0.06
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

def __sample_once__(src, nb, hm, cv, sz):
    """
    sample {sz} vertices centered at {cv} vertext in hemersphere {hm}
    across subject white matter stored in directory {src}, according
    to vertex neighborhood specified in {nb}.
    The white matter surfaces are supposedly in *wm.npz format.
    """
    ## mark vertices
    idx = [cv]
    mrk = set(idx)
    ## neighbor format: v_idx, n_nbr, nb[0], nb[1], ... nb[n_nbr-1]
    for i in xrange(sz):
        if len(idx) < sz:
            new = set(nb[idx[i]][2:]).difference(mrk)
            mrk.update(new)
            idx.extend(new)
        else:
            idx = idx[:sz]
            break

    sample = []
    for f, s in hlp.itr_fn(
            src, fmt = 'nc',
            flt = lambda w: w.endswith('wm.npz')):
        sample.append(np.load(f)[hm][idx])
    return np.array(sample)
    
def vtx_sample(src, dst = None, n = 7, sz = 10, seed = None):
    """ randomly sample regions from wm across
    all subjects in the source directory """
    import random
    random.seed(seed)

    if dst == None:
        dst = src
    if seed == None:
        seed = 0
    sfx = '{:02X}_{:04X}'.format(sz, seed)
    dst = pt.normpath(pt.join(dst, sfx))
    hlp.mk_dir(dst)

    ## get sample size, vertex count
    n = 2 ** n
    sz = 2 ** sz

    ## read vertex neighborhood information
    vnb = pt.join(src, __VNB)
    vnb = hlp.get_pk(vnb)

    ## choose hemersphere
    hem = [['lh', 'rh'][random.randint(0, 1)] for i in xrange(n)]

    ## choose center vertices
    cvx = {}
    cvx['lh'] = random.sample(xrange(len(vnb['lh'])), n)
    cvx['rh'] = random.sample(xrange(len(vnb['rh'])), n)
    cvx = [cvx[h][i] for i, h in zip(xrange(n), hem)]

    ## sample n sub surfaces
    for i in xrange(n):
        hm = hem[i]        # get hemesphere
        nb = vnb[hm]       # get neighborhood
        cv = cvx[i]        # get center vertex

        ## make file name, check overwiting
        fo = pt.join(dst, '{}{:05X}.npy'.format(hm, cv))
        if pt.isfile(fo):
            print fo, ": exists"
            continue

        ## sample, then save
        sp = __sample_once__(src, nb, hm, cv, sz)
        np.save(fo, sp)
        print fo, ": created"

def test():
    write_align_script(dst = '../tmp/align_vtx')
    vtx_asc2npz('../tmp/align_vtx', dst = '../tmp/wm_asc2npz')
    vtx_sample('../tmp/wm_asc2npz', '../tmp/wm_sample', seed = 120)
    
if __name__ == "__main__":
    pass
