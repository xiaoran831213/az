import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip
from shutil import copy as cp

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
    
def write_align_script(dst, ids = None, nodes = 4, psz = 8, hpc = True):
    """
    ids: list of subject to call align_value.sh, None mean use all subjects
    in $SUBJECTS_DIR
    dst: output destination of vertex alignment
    nodes: number of nodes to request from HPCC per submission.
    psz: process size - how much command to run on one nodes.
    """
    pfx = 'script'       # where to write scripts?
    hlp.mk_dir(pt.join(dst, pfx))
    cp('align_vtx.sh', pt.join(dst, pfx))
    dst=hlp.resolve_path(dst)

    ## read subject ids
    ids=__resolve__sid__(ids)

    ## container to hold the combined recon-all command
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    mem = nodes * 1.00            # for hpc
    wtm = psz * 0.08            # for hpc
    nbat = 0
    bsz = nodes * psz             # batch size
    cmd = pfx + '/align_vtx.sh -s {sid} -d . &>> {sid}.log\n'
    for i, s in enumerate(ids):
        j = i % bsz        # within batch index
        if j == 0:         # new batch
            f = open(pt.join(dst, pfx, fbat.format(nbat)), 'wb')
            if hpc:
                hlp.write_hpcc_header(
                    f, mem = mem, wtm = wtm, nodes = nodes)
                f.write('\n')
                
        k = j % psz        # within nodes index
        if k == 0:         # new nodes
            f.write('## node {:02d}\n'.format(j / psz))
            f.write('(\n')
            
        ## write new command
        f.write(cmd.format(sid=s))

        ## end of one nodes line
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
        bat = pt.join(pfx, fbat.format(i))
        f.write(frun.format(bat))
    f.close()
    hlp.chmod_x(f)

## type definitions
__NPY = np.dtype([
    ("x", '<f4'),
    ("y", '<f4'),
    ("z", '<f4'),
    ("spc", "<f4"),
    ("crv", "<f4"),
    ("slc", "<f4"),
    ("tck", "<f4")])
    
__ASC = [
    (0, "x", float),
    (1, "y", float),
    (2, "y", float),
    (3, "spc", float),
    (4, "crv", float),
    (5, "slc", float),
    (6, "tck", float)]

def __extract_neighbor_table__(dst = None, ovr = False):
    """
    extract vertex neighborhood table from the white matter surface
    of freesurfer atlas, and save them into python pickle format.
    """
    from subprocess import call
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

            ## create numpy array
            rs = np.array(rs, dtype = __NPY)

        ## the white matter surface in numpy
        np.savez_compressed(fo, lh=ls, rh=rs)
        print fo, ": created"
        
    ## extract vertex neighborhood table
    __extract_neighbor_table__(dst, ovr)

def test():
    write_align_script(dst = '../tmp/align_vtx')
    wm_asc2npz('../tmp/align_vtx', dst = '../tmp/wm_asc2npz')
    
if __name__ == "__main__":
    sid = __resolve__sid__()
    sdr = pt.expandvars('$SUBJECTS_DIR')

    for f, s in hlp.hpcc_iter(
        sid, '../hpc/align_vtx', npb = 4, qsz = 6, tpp = 0.20,
        mds = ['FreeSurfer/5.3.0'], lnk = ['align_vtx.sh'],
        pfx = ['export SUBJECTS_DIR={}'.format(sdr)],
            debug = False):

        f.write('./align_vtx.sh -s {sbj} > {sbj}.log\n'.format(sbj=s))
    
    pass
