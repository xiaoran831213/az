## comman helper functions
import cPickle
from glob import glob as gg
import os.path as pt
import os
import sys
import pdb
import shutil

def itr_fn(src = "", fmt = 'n', flt = None, drop = True):

    """
    filename iterator

    drop: drop list structure if only one file attribute
    is returned.
    format code:
    n: file name, N: absolute file name
    c: core name, C: absolute core name
    b: base name, B: absolute base name
    d: directory, D: absolute directory
    e: extension, E: absolute extension
    """
    src = resolve_path(src)
    if pt.isdir(src):
        src = pt.join(src, "*")

    if flt == None:
        flt = lambda w: True

    i = 0
    for fn in gg(src):
        if not flt(fn):
            continue
        rt = []
        for c in fmt:
            if c == 'i':
                r = i
            elif c == 'n':
                r = fn
            elif c == 'N':                # absolute filename
                r = pt.abspath(fn)
            elif c == 'C':                # absolute corename
                r = pt.abspath(fn).split('.')[0]
            elif c == 'c':                # ralative corename
                r = pt.basename(fn).split('.')[0]
            elif c == 'B':                # basename.extension
                r = pt.basename(pt.abspath(fn))
            elif c == 'b':                # basename
                r = pt.basename(fn)       
            elif c == 'D':                # absolute directory
                r = pt.dirname(pt.abspath(fn))
            elif c == 'd':                # relative directory
                r = pt.dirname(fn)
            elif c == 'e':                # extension(s)
                r = pt.basename(fn).split('.')[1:]
                if len(r) == 1:
                    r = r[0]
                if len(r) == 0:
                    r = None
            elif c == 'E':
                r = pt.basename(fn).split('.')[1:]
                if len(r) > 0:
                    r = r[-1]
                if len(r) == 0:
                    r = None
            else:
                continue
            rt.append(r)
        i += 1
        if drop and len(rt) == 1:
            yield rt[0]
        else:
            yield rt

def save_pgz(fo, s):
    """ save python object to gziped pickle """
    import gzip
    import cPickle
    with gzip.open(fo, 'wb') as gz:
        cPickle.dump(s, gz, cPickle.HIGHEST_PROTOCOL)

def load_pgz(fi):
    """ load python object from gziped pickle """
    import gzip
    import cPickle
    with gzip.open(fi, 'rb') as gz:
        return cPickle.load(gz)

def get_pk(src, idx = 0):
    """ get data from pickle """
    if pt.isdir(src):
        fn = gg(pt.join(src, "*"))[idx]
    else:
        fn = src

    with open(fn, 'rb') as fp:
        obj = cPickle.load(fp)

    print fn + ": fetched"
    return obj

def set_pk(obj, dst):
    if pt.isdir(dst):
        mk_dir(dst)
        fn = pt.join(dst, obj.__name__)
    else:
        mk_dir(pt.dirname(dst))
        fn = dst

    with open(dst, 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)
    print fn, ": dumpped"

def num_pk(src):
    return len(gg(src))

def itr_pk(src, fmt = ''):
    """ pickle file iterator """
    fmt = 'n' + fmt
    for rt in itr_fn(src, fmt, drop = False):
        with open(rt[0], 'rb') as pk:
            rt[0] = cPickle.load(pk)
        yield rt

def mk_dir(d):
    """ make deep folder """
    try:
        os.makedirs(d)
    except OSError as e:
        if not e.args[1] == 'File exists':
            raise e

def AUC(x, z):
    from sklearn.metrics import roc_auc_score
    import numpy
    x = x.reshape(x.shape[0], -1)
    z = z.reshape(z.shape[0], -1)
    s = numpy.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s.mean()

def resolve_path(path, full = False):
    path = pt.expanduser(path)
    path = pt.expandvars(path)
    if full:
        path = pt.abspath(path)
    return path

class hpcc_iter:
    def __init__(
            self,
            src, dst = None, npb = 3, ppn = 1, qsz = 1, mpn = 1, tpp = 0.1,
            mds = None, cpy = None, lnk = None, pfx = None, sfx = None, debug = False):
        """
        src: an iterable to supply commnand specific information
        dst: an directory to store generated HPCC scripts

        npb: nodes per batch
        qsz: queue size, number of processes to be lined up at each node
        ppn: processor per node. multiple processors per node will be helpful
        if the command to be wrapped has built in parallelism (e.g. OMP).
        mpn: memory per node, in GB.
        tpp: time per process, in Hour.

        mds: a list of module(s) to be loaded
        cpy: resources to be copied to the script environment (e.g. source coded).
        lnk: resources to be linked from the script environment (e.g. data files).
        pfx: batch prefix, a list of strings to write to the header of a batch.
        sfx: batch suffix, a list of strings to write to the footer of a batch.
        """
        self.itr = iter(src)
        
        self.dst = '/tmp/hpc' if dst is None else pt.normpath(dst)
        self.sdr = 'pbs'                            # script directory
        mk_dir('{}/{}'.format(self.dst, self.sdr))

        self.npb = 1 if npb is None else npb
        self.qsz = qsz
        self.ppn = 1 if ppn is None else ppn
        self.mpn = 1 if mpn is None else mpn
        self.tpp = 1 if tpp is None else tpp

        self.bsz = self.npb * self.qsz          # batch size

        self.mds = [] if mds is None else mds     # modules to load
        self.pfx = [] if pfx is None else pfx     # batch prefix
        self.sfx = [] if sfx is None else sfx     # batch surfix

        cpy = [] if cpy is None else cpy
        for r in cpy:
            shutil.copy(r, self.dst)

        lnk = [] if lnk is None else lnk
        for r in lnk:
            r = pt.abspath(pt.expandvars(pt.expanduser(r)))
            d = '{}/{}'.format(self.dst, pt.basename(pt.normpath(r)))
            if pt.exists(d):
                os.remove(d)
            os.symlink(pt.abspath(r),  d)

        ## iteration contex
        self.fo = open(os.devnull, 'wb')         # current script file
        self.np = 0                              # number of processes

        self.debug = debug

    def __write_hpcc_header__(self):
        ## batch index
        ibt = self.np / self.bsz
        
        ## open a new script file
        if self.debug:
            self.fo = sys.stdout
        else:
            fbt = '{}/{}/{:03d}.sh'.format(self.dst, self.sdr, ibt)
            dst = self.dst
            dst = pt.expanduser(dst)
            dst = pt.expandvars(dst)
            self.fo = open(fbt, 'wb')

        self.fo.write('#!/bin/bash -login\n')

        ## -------- resource requsition -------- ##
        ## nodes and processor per node
        self.fo.write('#PBS -l nodes={}:ppn={}\n'.format(self.npb, self.ppn))
        ## walltime
        wtm = self.qsz * self.tpp
        hh = int(wtm);
        wtm -= hh
        mm = int(wtm * 60);
        self.fo.write('#PBS -l walltime={:02d}:{:02d}:00\n'.format(hh, mm))
        ## memory
        mem = self.npb * self.mpn
        self.fo.write('#PBS -l mem={}M\n'.format(int(mem*1024)))
        ## others
        self.fo.write('#PBS -j oe\n')       # combin stdout and stderr
        self.fo.write('#PBS -V\n')          # copy environment vars
        self.fo.write('#PBS -o std\n')      # store captured stdout

        self.fo.write('\n')

        ## module loading
        if len(self.mds) > 0:
            fmdl = 'module load {}\n'
            self.fo.writelines([fmdl.format(m) for m in self.mds])
            self.fo.write('\n')

        ## working directory
        self.fo.write('[ -n "$PBS_O_WORKDIR" ] && cd "$PBS_O_WORKDIR"\n\n')

        ## batch prefix
        if len(self.pfx) > 0:
            self.fo.write('\n'.join(self.pfx))
            self.fo.write('\n\n')

        ## batch index
        self.fo.write(
            '## -------- batch {:03d} -------- ##\n'.format(ibt))

    def __write_hpcc_footer__(self):
        ## wait sub processes running on multiple nodes to finish
        if self.npb > 1:
            self.fo.write('wait')
        self.fo.write('\n')

        ## batch suffix
        if len(self.sfx) > 0:
            self.fo.write('\n'.join(self.sfx))
            self.fo.write('\n')

        ## close the script file
        if self.fo is not sys.stdout:
            self.fo.close()

    def __write_node_header__(self):
        if self.npb > 1:
            i_node = (self.np % self.bsz) / self.qsz
            self.fo.write('## node {:02d}\n'.format(i_node))
            self.fo.write('(\n')

    def __write_node_footer__(self):
        if self.npb > 1:
            self.fo.write(')&\n\n')

    def __write_submiter__(self):
        nbt = self.np / self.bsz
        if self.np % self.bsz > 0:
            nbt += 1

        fsb = '({} {}/{{:03d}}.sh)&\n'.format('qsub', self.sdr)
        if self.debug:
            f = sys.stdout
        else:
            f = open('{}/sub.sh'.format(self.dst), 'wb')

        f.write('#!/bin/bash\n')
        f.write('cd "`dirname $0`"\n')
        f.write('test -d std || mkdir std\n')
        for i in xrange(nbt):
            f.write(fsb.format(i))
            if (i + 1) % 8 == 0:
                f.write('wait\n')
        f.write('wait\n')

        if f is not sys.stdout:
            f.close()
            chmod_x(f)

    def __iter__(self):
        return self

    def next(self):
        try:
            cmdobj = self.itr.next()       # command specific object
        except StopIteration:
            self.__write_node_footer__()
            self.__write_hpcc_footer__()
            self.__write_submiter__()
            raise StopIteration()

        if self.np % self.bsz is 0:        # new batch
            self.__write_node_footer__()   # wrap up previous script  
            self.__write_hpcc_footer__()   
            self.__write_hpcc_header__()
            self.__write_node_header__()
        elif self.np % self.qsz is 0:      # new node
            self.__write_node_footer__()   # wrap up previous node
            self.__write_node_header__() 

        self.np += 1
        return self.fo, cmdobj

def chmod_x(fi):
    ## make the script executable
    import stat as st
    fn = fi.name if isinstance(fi, file) else fi
    mode = os.stat(fn).st_mode
    os.chmod(fn, mode|st.S_IXUSR|st.S_IXGRP|st.S_IXOTH)
    
if __name__ == "__main__":
    hi = hpcc_iter(
        range(20), dst = '../tmp', qsz = 3, npb = 1,
        mds = ['R/3.1.0'], lnk = ['align_vtx.sh', '../dat'], debug = True)
    for fo, cm in hi:
        print fo, cm
    pass
