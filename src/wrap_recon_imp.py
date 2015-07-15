import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import sys
import hlp
import pdb

__NPY = np.dtype([
    ("img", "<i4"),
    ("sbj", "a16"),
    ("grp", "a16"),
    ("sex", "a08"),
    ("age", "<i1"),
    ("vst", "<i1"),
    ("mod", "a16"),
    ("dsc", "a32"),
    ("typ", "a16"),
    ("acq", "a16"),
    ("fmt", "a08")])
    
__CSV = [
    ("Image Data ID", int),
    ("Subject", str),
    ("Group", str),
    ("Sex", str),
    ("Age", int),
    ("Visit", int),
    ("Modality", str),
    ("Description", str),
    ("Type", str),
    ("Acq Date", str),
    ("Format", str)]

def __read_manifest__(manifest_csv):
    """ read ADNI downloaded image manifest_csv CSV file"""
    fi = open(manifest_csv, 'rb')
        
    ## prepare csv reader, deduce field name from header
    ls = []
    for line in csv.DictReader(fi,fieldnames = None):
        ls.append(tuple([t(line[k]) for k, t in __CSV]))
    mf = np.array(ls, dtype = __NPY)
    fi.close()
    return mf

def write_recon_imp(csv, dst, nodes = 4, mpn = 1, ppn = 1, qsz = 4, tpp = 0.10, hpc = 0):
    """
    write HPCC script for freesurfer image imporation.
    """
    pfx = 'script'
    hlp.mk_dir(pt.join(dst, pfx))
        
    ## load manifest CSV and sort by subject id
    mf = __read_manifest__(csv)
    mf = mf[mf['sbj'].argsort()]
    
    ## iterate MRI scans, group them into
    from itertools import groupby
    
    ## write hpcc script now
    nbat = 0
    fbat = '{:03d}.qs' if hpc else '{:03d}.sh'
    fimp = '-i $ADNI_IMG/{img:08d}.nii'
    fcmd = 'recon-all {imp} -s {sbj} &> {sbj}.log\n'
    nbat = 0
    bsz = nodes * qsz             # batch size
    f = None                      # script file
    for i, sbj_grp in enumerate(groupby(mf, lambda w: w['sbj'])):
        sbj, grp = sbj_grp[0], sbj_grp[1]
        j = i % bsz               # within batch index
        if j == 0:                # new batch
            f = open(pt.join(dst, pfx, fbat.format(nbat)), 'wb')
            hlp.write_hpcc_header(
                f,
                mem = nodes * mpn,
                wtm = qsz * tpp,
                nodes = nodes)
            f.write('\n')

        ## new node 
        if j % qsz is 0 and nodes > 1:
            f.write('## node {:02d}\n'.format(j/qsz))
            f.write('(\n')

        ## go through all MRI scan of one subject to write FreeSurfer
        ## image importing command
        ## imp = " ".join([fimp.format(img=mri['img'].item()) for mri in grp])
        imp = fimp.format(img = iter(grp).next()['img'].item())
        f.write(fcmd.format(imp=imp, sbj=sbj))

        ## end of the node 
        if (j + 1) % qsz is 0 and nodes > 0:
            f.write(')&\n\n')
        
        ## end of one batch
        if (i + 1) % bsz == 0:
            if nodes > 1:
                f.write('wait\n')
            nbat += 1
            f.close()

    # the left over
    if f is not None and not f.closed:
        if nodes > 1:
            if (j + 1) % qsz is not 0:
                f.write(')&\n\n')
            f.write('wait\n')
        nbat += 1
        f.close()

    ## write submition script
    f = open(pt.join(dst, 'sub.sh'), 'wb')
    f.write('#!/bin/bash\n')
    f.write('cd "`dirname $0`"\n')
    ipt = 'qsub' if hpc else 'sh'         # choose interpreter
    for i in xrange(nbat):
        f.write('{} {}/{:03d}.sh\n'.format(ipt, pfx, i))
    f.close()
    hlp.chmod_x(f)

def main():
    write_recon_imp('../dat/ADNI_G800.csv', '../hpc/recon_imp')
        
if __name__ == "__main__":
    pass
