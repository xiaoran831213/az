import os
import os.path as pt
import csv
import mri.work as mw
import numpy as np
import hlp
import pdb
import rdm.sda as sda


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
        use = 'dat/img/' + str(lb)
        mw.vtx2vox(npy, vox, ovr = 0, flt = lambda v: v['lbl'] == lb)
        mw.sfr2vlm(vox, vlm, ovr = 0, dim = 48)
        mw.vlm2vmk(vlm, vmk, ovr = 0)
        mw.pack(vmk, use, ovr = 0)

def train_sda(src, dst, ovr = 0):
    hlp.mk_dir(dst)
    for x, lb in hlp.itr_pk('dat/img', fmt = 'b'):
        fo = pt.join('dat/sda', lb)
        if pt.exists(fo) and not ovr:
            print fo, "exists"
            continue

        x = x.reshape(x.shape[0], -1)
        s = sda.SDA(n_vis = x.size/x.shape[0], n_hid = (270, 90, 30, 10))

        print "pre-train:"
        sda.pre_train(s, x)
        z = s.f_pred()(x)
        print "AUC after pre-train: ", hlp.AUC(x, z)
        print hlp.AUC(x, z)

        print "fine-tune:"
        sda.fine_tune(s, x)
        z = s.f_pred()(x)
        print "AUC after fine-tune: ", hlp.AUC(x, z)
        hlp.set_pk(s, fo)

def encode_img(src, dst, ovr = 0):
    hlp.mk_dir(dst)
    for x, lb in hlp.itr_pk(src, fmt = 'b'):
        fo = pt.join(dst, lb)
        if pt.exists(fo) and not ovr:
            print fo, "exists"
            continue

        x = x.reshape(x.shape[0], -1)
        s = hlp.get_pk(pt.join('dat/sda', lb))

        y = s.f_encode()(x)
        np.savetxt(fo, y)
        print fo, "created"

def ftmp():
    ssn = mw.pack('dat/vmk/1003', 'tmp/pck/1003', ovr = 1)
    with open('tmp/ssn.1003', 'w') as f:
        for sn in ssn:
            f.write(sn + '\n')

    ssn = mw.pack('dat/vmk/1035', 'tmp/pck/1035', ovr = 1)
    with open('tmp/ssn.1035', 'w') as f:
        for sn in ssn:
            f.write(sn + '\n')

    ssn = mw.pack('dat/vmk/2003', 'tmp/pck/2003', ovr = 1)
    with open('tmp/ssn.2003', 'w') as f:
        for sn in ssn:
            f.write(sn + '\n')

    ssn = mw.pack('dat/vmk/2035', 'tmp/pck/2035', ovr = 1)
    with open('tmp/ssn.2035', 'w') as f:
        for sn in ssn:
            f.write(sn + '\n')

if __name__ == "__main__":
    pass

