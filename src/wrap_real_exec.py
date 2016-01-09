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
    hlp.mk_dir(dst)

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
    ## wms = sorted(gg(pt.join(img, '*.rds')))
    wms = ['lh31']
    chs = reduce(lambda x,y:x+y, ((1+c, 24-c) for c in xrange(0, 12)))
    for w, c in [(w, c) for w in wms for c in chs]:
        wm = pt.splitext(pt.basename(w))[0]
        ch = 'ch{:02d}'.format(c)
        fo = pt.join(dst, '{}_{}.rds'.format(wm, ch))
        fl = pt.join(dst, '{}_{}.log'.format(wm, ch))
        if pt.isfile(fo) and not ovr:
            print fo, ': exists'
        elif pt.isfile(fl) and not ovr:
            print fl, ': exists'
        else:
            items.append((wm, ch))

    ## write commands
    cmd = 'time ./run_mix.R img gno -w {w} -c {c} -d {w}_{c}.rds &> {w}_{c}.log\n'
    for fo, it in hlp.hpcc_iter(
            items, dst, npb=1, ppn= 4, mpn=4, tpp=3.0, qsz = 1,
            mds=['R/3.1.0'], lnk=['.', '../dat'],
            cpy=['run_mix.R'],
            debug=0):

        ## write command for one processor
        wm, ch = it[0], it[1]
        fo.write(cmd.format(w=wm, c=ch))
        print 'write:', wm, ch

def test():
    write_real_exec('$AZ_AENC', '$AZ_WGS.bin', dst = '$AZ_REAL', ovr = 0)

if __name__ == '__main__':
    import os
    import sys
    pass
