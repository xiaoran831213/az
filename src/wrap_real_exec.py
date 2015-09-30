import pdb
import os
import os.path as pt
import hlp

def write_real_exec(img, gno, dst = None, ovr = 0):
    """
    prepare data for real analysis
    """
    ## resove source and target directories
    img = pt.expandvars(pt.expanduser(img))
    gno = pt.expandvars(pt.expanduser(gno))
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

    from glob import glob as gg
    items = []
    for i, c in [(i, c) for i in gg(pt.join(img, '*.rds')) for c in xrange(1, 25)]:
        im = pt.splitext(pt.basename(i))[0]
        ch = '{:02d}'.format(c)
        fo = pt.join(dst, '{}_{}.rds'.format(im, ch))
        fl = pt.join(dst, '{}_{}.log'.format(im, ch))
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        elif pt.isfile(fl) and not ovr:
            print fl, ': exists'
        else:
            items.append((im, ch))
    items.sort()
    pdb.set_trace()
    ## write commands
    cmd = './run_mix.R img gno -i {i} -c {c} -d {i}_{c}.rds &> {i}_{c}.log\n'
    for fo, it in rut_hlp.hpcc_iter(
            items, dst, npb=1, ppn= 4, mpn=2, tpp=3.0, qsz = 1,
            mds=['R/3.1.0'], lnk=['.'], 
            debug=0):

        ## write command for one processor
        im, ch = it[0], it[1]
        fo.write(cmd.format(i=im, c=ch))

def test():
    write_real_exec('$AZ_AENC', '$AZ_WGS.bin', dst = '.', ovr = 0)

if __name__ == '__main__':
    ## add project root to python path
    import os
    import sys
    pass
