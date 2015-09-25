import pdb
import numpy as np
import os
import os.path as pt

def write_pre_train(src, dst, ovr = 0):
    """
    train SDA with WM samples in {src}
    zsd: region size and random seed used to sample WM, which decides the folder name
    of the samples, and the default destination folder for SDA
    """
    src = pt.expandvars(pt.expanduser(src))
    dst = pt.expandvars(pt.expanduser(dst))
    rut_hlp.mk_dir(dst)

    ## gather WM surface samples, also check existing output
    sfs = []
    flt = lambda w: w.endswith('npz') or w.endswith('pgz')
    for sf in rut_hlp.itr_fn(src, 'c', flt):
        fo = pt.join(dst, sf + '.pgz')
        fl = pt.join(dst, sf + '.log')
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        elif pt.isfile(fl) and not ovr:
            print fl, ': exists'
        else:
            sfs.append(sf)

    ## write commands
    cmd = 'time python wm_pre_train.py {t} {s} . &>{t}.log\n'
    for fo, sf in rut_hlp.hpcc_iter(
            sfs, dst, npb=1, ppn= 4, mpn=4, tpp=4, qsz = 1,
            mds=['NumPy', 'R/3.1.0'],
            lnk=['rdm', 'hlp.py'],
            cpy=['wm_pre_train.py'],
            pfx=['export MKL_NUM_THREADS={}'.format(4)],
            debug=False):

        ## write command for one processor
        fo.write(cmd.format(s=src, t=sf))

def test():
    write_pre_train('$AZ_APEC', '$AZ_APTN')
    pass

if __name__ == '__main__':
    ## add project root to python path
    import os
    import sys
    pass
