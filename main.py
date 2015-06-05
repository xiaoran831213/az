import csv
import os
import mri.work as mw
import numpy as np
import hlp
import pdb

def test():
    import mri.work as nw
    n0 = []
    for sf, sn in hlp.itr_pk('dat/vox/2035/', bsn = True):
        #print sn, sf.shape[0]
        n0.append(sf.shape[0])

    n1 = []
    for vl, sn in hlp.itr_pk('dat/vlm/2035/', bsn = True):
        #print sn, np.count_nonzero(vl['lbl'])
        n1.append(np.count_nonzero(vl['lbl']))

    n0 = np.array(n0)
    n1 = np.array(n1)
    lost = (n0 - n1)/np.asarray(n0, dtype='<f4')
    return n0, n1, lost

def proc_image_data():
    # from mri.tar2csv import tar2csv as tar2csv
    # tar2csv('raw/vtx.tar.gz', 'dat/csv')
    labels = np.unique(hlp.get_pk('dat/npy')['lbl']).tolist()
    labels.remove(-1)
    for lb in labels:
        npy = 'dat/npy'
        vox = 'dat/vox/' + str(lb)
        vlm = 'dat/vlm/' + str(lb)
        vmk = 'dat/vmk/' + str(lb)
        use = 'dat/use/' + str(lb)
        mw.vtx2vox(npy, vox, ovr = 0, flt = lambda v: v['lbl'] == lb)
        mw.sfr2vlm(vox, vlm, ovr = 0, dim = 48)
        mw.vlm2vmk(vlm, vmk, ovr = 0)
        mw.pack(vmk, use, ovr = 0)

if __name__ == "__main__":
    reload(mw)
    pass
