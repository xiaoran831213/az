import csv
import os
import mri.work as mw
import numpy as np
import hlp

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
    mw.vtx2vox('dat/npy', 'dat/vox/1003', ovr = 0, flt = lambda v: v['lbl'] == 1003)
    mw.sfr2vlm('dat/vox/1003', 'dat/vlm/1003', ovr = 1, dim = 48)
    mw.vlm2vmk('dat/vlm/1003', 'dat/vmk/1003', ovr = 1)
    mw.pack('dat/vmk/1003', 'dat/use/1003', ovr = 1)

    mw.vtx2vox('dat/npy', 'dat/vox/1035', ovr = 0, flt = lambda v: v['lbl'] == 1035)
    mw.sfr2vlm('dat/vox/1035', 'dat/vlm/1035', ovr = 1, dim = 48)
    mw.vlm2vmk('dat/vlm/1035', 'dat/vmk/1035', ovr = 1)
    mw.pack('dat/vmk/1035', 'dat/use/1035', ovr = 1)
    
    mw.vtx2vox('dat/npy', 'dat/vox/2003', ovr = 0, flt = lambda v: v['lbl'] == 2003)
    mw.sfr2vlm('dat/vox/2003', 'dat/vlm/2003', ovr = 1, dim = 48)
    mw.vlm2vmk('dat/vlm/2003', 'dat/vmk/2003', ovr = 1)
    mw.pack('dat/vmk/2003', 'dat/use/2003', ovr = 1)
    
    mw.vtx2vox('dat/npy', 'dat/vox/2035', ovr = 0, flt = lambda v: v['lbl'] == 2035)
    mw.sfr2vlm('dat/vox/2035', 'dat/vlm/2035', ovr = 1, dim = 48)
    mw.vlm2vmk('dat/vlm/2035', 'dat/vmk/2035', ovr = 1)
    mw.pack('dat/vmk/2035', 'dat/use/2035', ovr = 1)
    
if __name__ == "__main__":
    reload(mw)
    pass
