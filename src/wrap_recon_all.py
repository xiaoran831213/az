import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import sys
import hlp
# -subjid <subjid>
# -<process directive>

# Fully-Automated Directive:
# -all           : performs all stages of cortical reconstruction
# -autorecon-all : same as -all

##  -noappend     : start new log and status files instead of appending
##  -no-isrunning : do not check whether this subject is currently being processed
def write_recon_all(dst, nodes = 1, ppn = 1, qsz = 1, mpn = 4, tpp = 24, hpc = 1):
    """
    write HPCC PBS script to call recon-all on imported subjects
    dst: destination directory to store command output and hpcc scripts
    nodes: the number of nodes to request from HPCC for one batch
    ppn: processors per node, useful when the command has built in parallelism
    qsz: queue size - number of processes lining up on a node
    mpn: memory per node (in GB)
    tpp: time per process (in Hour)
    """
    pfx = 'script'
    hlp.mk_dir(pt.join(dst, pfx))

    ## the source directory is freesurfer's subject directory
    src = hlp.resolve_path('$SUBJECTS_DIR', full = True)

    ## container to hold the combined recon-all command
    fcmd = 'recon-all -all -s {0} -noappend -no-isrunning &> {0}.log\n'
    fbat = '{dst}/{pfx}/{bat:03d}.sh'
    bsz = nodes * qsz
    nbat = 0
    flt = lambda fn: pt.isdir(fn)
    for i, s in enumerate(hlp.itr_fn(src, fmt = 'b', flt = flt)):
        ## i: subject index
        ## s: subject ID
        j = i % bsz        # within batch index
        if j is 0:         # new batch
            f = open(fbat.format(dst=dst, pfx=pfx, bat=nbat), 'wb')
            hlp.write_hpcc_header(
                f,
                mem = nodes * mpn,
                wtm = tpp * qsz,
                nodes = nodes)
            f.write('\n')

        ## a new node
        k = j % nodes                     # within node index
        if k is 0 and nodes > 1:
            f.write('## node {:02d}\n'.format(j / qsz))
            f.write('(\n')
            
        ## write new command
        f.write(fcmd.format(s))

        ## end of the node 
        if (j + 1) % qsz is 0 and nodes > 1:
            f.write(')&\n\n')
        
        ## end of one batch
        if (i + 1) % bsz == 0:
            if nodes > 1:
                f.write('wait\n')
            nbat += 1
            f.close()

    # the left over
    if f in locals() and not f.closed:
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
    write_recon_all(dst = '../hpc/recon_all')

if __name__ == "__main__":
    pass
