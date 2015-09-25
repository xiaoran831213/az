import pdb
import os
import os.path as pt

def write_simu_exec(img, gno, dst = None, ovr = 0):
    """
    prepare data for simulation studies.
    gaussian blur the source WM surfaces.
    """
    ## resove source and target directories
    img = pt.expandvars(pt.expanduser(img))
    gno = pt.expandvars(pt.expanduser(gno))
    if dst is None:
        dst = '/tmp/sim'
    else:
        dst = pt.expandvars(pt.expanduser(dst))
    rut_hlp.mk_dir(dst)

    ## use the symbolic link to the source folder to shorten the script
    lnk = pt.join(dst, 'img')
    if pt.exists(lnk):
        os.remove(lnk)
    os.symlink(img, lnk)

    lnk = pt.join(dst, 'gno')
    if pt.exists(lnk):
        os.remove(lnk)
    os.symlink(gno, lnk)

    itr = []
    for it in xrange(500):
        it_str = '{:03d}'.format(it)
        fo = pt.join(dst, it_str + '.rds')
        fl = pt.join(dst, it_str + '.log')
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        elif pt.isfile(fl) and not ovr:
            print fl, ': exists'
        else:
            itr.append(it)

    ## write commands
    cmd = './sim_mix.R img gno -i 20 -d {i:03d} -v tck &> {i:03d}.log\n'
    for fo, it in rut_hlp.hpcc_iter(
            itr, dst, npb=1, ppn= 4, mpn=2, tpp=3.0, qsz = 1,
            mds=['R/3.1.0'], lnk=['.'], 
            debug=0):

        ## write command for one processor
        fo.write(cmd.format(i=it))

def test():
    write_simu_exec('$AZ_SM2', '$AZ_WGS.bin', dst = '../sim/N1_5', ovr = 0)

if __name__ == '__main__':
    ## add project root to python path
    import os
    import sys
    if not os.environ['AZ_PRJ'] in sys.path:
        sys.path.insert(0, os.environ['AZ_PRJ'])
    import src.hlp as rut_hlp
    pass
