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
def write_recon_all_script(
    src = None, flt = None, dst = None, ncpu = 4, ncmd = 8):
    """
    write HPCC PBS script to call recon-all on imported subjects
    """
    if dst == None:
        dst = "hpc/recon-all"
    hlp.mk_dir(dst)

    if src == None:
        src = '$SUBJECTS_DIR'
    src = hlp.resolve_path(src, full = True)

    if flt == None:
        flt = lambda fn: pt.isdir(fn)

    ## container to hold the combined recon-all command
    cmd = 'recon-all -all -s {0} -noappend -no-isrunning 1>{0}.out 2>{0}.err || echo -n\n'
    bat = '{}/{:03d}.qs'
    mem = ncpu * 4.0
    batch_size = ncpu * ncmd
    wtm = 24
    nbat = 0
    for i, s in enumerate(hlp.itr_fn(src, fmt = 'b', flt = flt)):
        j = i % batch_size        # within batch index
        if j == 0:                # new batch
            f = open(bat.format(dst, nbat), 'wb')
            hlp.write_hpcc_header(
                f, mem = mem, walltime = wtm, nodes = ncpu)
            f.write('\n')
            
        k = j % ncmd              # within cpu index
        if k == 0:                # new cpu
            f.write('## node {:02d}\n'.format(j / ncmd))
            f.write('(\n')
            
        ## write new command
        f.write(cmd.format(s))

        ## end of one cpu line
        if (j + 1) % ncmd == 0:
            f.write(')&\n\n')
        
        ## end of one batch
        if (i + 1) % batch_size == 0:
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
    for ibat in xrange(nbat):
        f_bat = pt.abspath(bat.format(dst, ibat))
        f.write('qsub {}\n'.format(f_bat))
        
def main():
    import os
    import sys
    if len(sys.argv) < 3:
        print "Usage: ", sys.argv[0], "<sbj_dir> <*dst>" 
        return
    
    src = sys.argv[1]
    dst = None
    if len(sys.argv) > 2:
        dst = sys.argv[2]

    write_recon_all_script(src, dst = None)

def test():
    write_recon_all_script(
        src = None, dst = '../tmp', ncpu = 1, ncmd = 1)
        
if __name__ == "__main__":
    pass
