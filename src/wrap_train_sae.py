import pdb
import numpy as np
import os
import os.path as pt

def write_train_sae(src, dst = None, ovr = 0):
    """
    train SDA with WM samples in {src}
    zsd: region size and random seed used to sample WM, which decides the folder name
    of the samples, and the default destination folder for SDA
    """
    ## resove source and target directories
    src = pt.expandvars(pt.expanduser(src))
    if dst is None:
        dst = pt.join(pt.dirname(src), 'trained_sae')
    else:
        dst = pt.expandvars(pt.expanduser(dst))
    rut_hlp.mk_dir(pt.join(dst))

    ## gather WM surface samples, also check existing output
    sfs = []
    for sf in rut_hlp.itr_fn(src, 'c', lambda w: w.endswith('npz')):
        fo = pt.join(dst, sf + '.pgz')
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        else:
            sfs.append(sf)

    ## write commands
    tsk = 'tsk/{t}.pk'
    cmd = 'time python rdm/train_sda.py {w} &>{t}.log\n'
    for fo, sf in rut_hlp.hpcc_iter(
            sfs, dst, npb=1, ppn= 8, mpn=.75, tpp=2, qsz = 2,
            mds=['NumPy', 'R/3.1.0'],
            lnk=['rdm'],
            pfx=['export MKL_NUM_THREADS={}'.format(8)],
            debug=False):

        ## save the working material specification for one processor
        whr = pt.join(dst, tsk).format(t=sf)
        ## src: wm sample directory, dst: distination directory
        ## wms: wm sample id
        wrk = {
            'src' : pt.abspath(src),            # source directory
            'dst' : pt.abspath(dst),            # target direcotry
            'wms' : sf}                         # white matter surface
        rut_hlp.set_pk(wrk, whr)
                
        ## write command for one processor
        fo.write(cmd.format(w=whr, t=sf))

def test():
    #write_train_sae('$AZ_SP1', dst = '$AZ_EC1', ovr = 1)
    write_train_sae('$AZ_SP1', dst = '$AZ_EC2', ovr = 0)
    pass

if __name__ == '__main__':
    ## add project root to python path
    import os
    import sys
    if not os.environ['AZ_PRJ'] in sys.path:
        sys.path.insert(0, os.environ['AZ_PRJ'])
    import src.hlp as rut_hlp
    pass
