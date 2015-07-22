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

    fo = pt.join(dst, 'vtx_nbr.ppk')
    if pt.isfile(fo) and not ovr:
        print fo, ': exists'
        return
    
    dc = {}
    s = pt.expandvars("$FREESURFER_HOME/subjects/fsaverage/surf")
    for h in ('lh', 'rh'):
        sfr = '{sf}/{hm}.white'.format(hm=h, sf=s)
        asc = '/tmp/{hm}.nbr.asc'.format(hm=h)
        cmd = ['mris_convert', '-v', sfr, asc]
        print cmd
        call(cmd)

        hem = []
        with open(asc) as f:
            for line in f:
                hem.append(tuple([int(w) for w in line.split()]))
        dc[h] = hem

    ## save the vertex neighbor table
    hlp.set_pk(dc, fo)


def __write_avg_sphere_cord__(dst = None):
    """
    write sphere coordination of fsaverage.
    both xyz and phe/theta are written.
    """
    from subprocess import call
    if dst == None:
        dst = "."
    else:
        hlp.mk_dir(dst)

    dt = np.dtype({
        'names':['x', 'y', 'z', 't', 'p'],
        'formats':['<f8', '<f8', '<f8', '<f8', '<f8']})
    sf = pt.expandvars("$FREESURFER_HOME/subjects/fsaverage/surf/{hm}.sphere.reg")
    fo = pt.join(dst, 'avg_sphere.npz')
    dc = {}
    for h in ('lh', 'rh'):
        sfr = sf.format(hm=h)
        tmp = '/tmp/{}.sphere.asc'.format(h)
        cmd = ['mris_convert', sfr, tmp]

        ## save ascii formated XYZ sphere coordinate to temporary files
        print cmd
        call(cmd)

        ## read into Numpy array
        f = open(tmp, 'rb')

        ## the 2nd line of the ascii table is vertex and triangle count,
        f.readline()                      # skip the first line
        nv, nt = [int(n) for n in f.readline().split(' ')]

        ## only vertex is need here, so skip_footer = nt
        ## the forth column is useless here
        s = np.genfromtxt(
            f, skip_footer=nt, names=['x', 'y', 'z', 'v'],
            usecols=['x', 'y', 'z'])
        f.close()

        ## calculate latitude and longitude coordinate
        x, y, z = [s[a] for a in 'xyz']         # xyz
        r = np.sqrt(x**2 + y**2 + z**2)         # radius

        s = np.ndarray(s.size, dtype = dt)
        s['x'] = x
        s['y'] = y
        s['z'] = z
        s['p'] = np.pi/2.0 - np.arccos(z/r)        # latitude, phi
        s['l'] = np.arctan2(y, x)                  # longitude, lambda
        dc[h] = s

    np.savez_compressed(dst, *dc)
    return dc

def main():
    sid = __resolve__sid__()
    sdr = pt.expandvars('$SUBJECTS_DIR')
    dst = '../hpc/align_vtx'

    for f, s in hlp.hpcc_iter(
            sid, dst, npb = 4, qsz = 6, tpp = 0.20,
            mds = ['FreeSurfer/5.3.0'],
            lnk = ['align_vtx.sh', 'wm_asc2npy.py'],
            pfx = ['export SUBJECTS_DIR={}'.format(sdr)],
            debug = False):
        f.write('./align_vtx.sh -s {sbj} > {sbj}.log\n'.format(sbj=s))
    __extract_neighbor_table__(dst, ovr = True)
    
if __name__ == "__main__":
    pass
