import pdb
import os
import os.path as pt

def write_simu_prep(src, dst = None, ovr = 0):
    """
    prepare data for simulation studies.
    gaussian blur the source WM surfaces.
    """
    ## resove source and target directories
    src = pt.expandvars(pt.expanduser(src))
    if dst is None:
        dst = pt.join(pt.dirname(src), 'trained_sae')
    else:
        dst = pt.expandvars(pt.expanduser(dst))
    rut_hlp.mk_dir(dst)

    ## use the symbolic link to the source folder to shorten the script
    src_lnk = pt.join(dst, 'wms')
    if pt.exists(src_lnk):
        os.remove(src_lnk)
    os.symlink(src, src_lnk)

    ## gather WM surface samples, also check existing output
    sfs = []
    for sf in rut_hlp.itr_fn(src, 'c', lambda w: w.endswith('rds')):
        fo = pt.join(dst, sf + '.rds')
        fl = pt.join(dst, sf + '.log')
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        elif pt.isfile(fl) and not ovr:
            print fl, ': exists'
        else:
            sfs.append(sf)

    ## write commands
    cmd = 'src/img.R wms {sn} &> {sn}.log\n'
    for fo, sf in rut_hlp.hpcc_iter(
            sfs, dst, npb=4, ppn= 4, mpn=1, tpp=.25, qsz = 4,
            mds=['R/3.1.0'], lnk=['.'], debug=0):

        ## write command for one processor
        fo.write(cmd.format(sn=sf))

def test():
    write_simu_prep('$AZ_EC2', dst = '$AZ_SM2', ovr = 0)
    pass

if __name__ == '__main__':
    ## add project root to python path
    import os
    import sys
    if not os.environ['AZ_PRJ'] in sys.path:
        sys.path.insert(0, os.environ['AZ_PRJ'])
    import src.hlp as rut_hlp
    pass
