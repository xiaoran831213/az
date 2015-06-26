import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import hlp
from itertools import izip

def __wm_neighbor__(nb, cv, sz):
    """
    According to vertex neighborhood {nb}, sample {sz}
    neighbors around the center vertex {cv}
    return vertex indices including the center one
    """
    ## mark vertices
    idx = [cv]
    mrk = set(idx)
    ## neighbor format: v_idx, n_nbr, nb[0], nb[1], ... nb[n_nbr-1]
    for i in xrange(sz):
        if len(idx) < sz:
            new = set(nb[idx[i]][2:]).difference(mrk)
            mrk.update(new)
            idx.extend(new)
        else:
            return idx[:sz]

def __sample_wm__(src, dst, hms, cvs):
    """
    given a list of center vertices {cvs}, and hemispheres {hms},
    pick regions from WM surface across subject in {src}
    """
    ## read vertex neighborhood information
    vnb = hlp.get_pk(pt.join(src, 'wm.vtxnbr.ppk'))
    nbs = (vnb[h] for h in hms)

    lvi, lfo = [], []  
    ## gather lists of neighbor indices and output file
    ## names for each center vertex
    for hm, nb, cv in izip(hms, nbs, cvs):
        ## make output file name, check overwiting
        fo = pt.join(dst, '{}{:05X}.npy'.format(hm, cv))
        if pt.isfile(fo):
            print fo, ": exists"
            continue
        lfo.append(fo)

        ## neighboring vertex indices
        lvi.append(__wm_neighbor__(nb, cv, sz))

    ## make lists of sampled surfaces and subject IDs
    lsf = [[]] * len(lvi)
    lsb = []
    print 'fetch %d WM areas from each subject:' % len(lvi) 
    for fn, sb in hlp.itr_fn(src, 'nc', flt = lambda w: w.endswith('wm.npz')):
        print sb
        wm = np.load(fn)
        for si, hm, vi in izip(xrange(len(lvi)), hms, lvi):
            lsf[si].append(wm[hm][vi])
        lsb.append(sb)
        
    for sf, fo in izip(lsf, lfo):
        np.save(fo, np.array(sf))
        print fo, ": created"

    with open(pt.join(dst, 'sbj.txt'), 'wb') as f:
        f.write("\n".join(lsb))

if __name__ == "__main__":
    import sys
    os.chdir('../tmp/wm_sample/0A_0078')

    if len(sys.argv) < 4:
        src = '/hd2/study/az/tmp/wm_asc2npz'
        dst = '.'
        h_c = 'script/004_03.pk'
    else:
        src = sys.argv[1]
        dst = sys.argv[2]
        h_c = sys.argv[3]
    hem, cvs = hlp.get_pk(h_c)
    pdb.set_trace()

    pass










