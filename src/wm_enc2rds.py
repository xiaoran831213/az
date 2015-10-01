import os.path as pt
import pdb
import numpy as np
import hlp

def __np2r__(x, dn):
    """ Numpy array to R array """
    import rpy2
    import rpy2.robjects as robjs
    import rpy2.rlike.container as rlc

    ## R only accept list based arguments
    shp = list(x.shape)
    typ = x.dtype.name

    ## transpose is needed since R grows the left most dimension fastest
    x = x.transpose().flatten().tolist()

    ## directly involking R.array without vector wrapping results in
    ## an array of lists of size 1
    if typ.startswith('int'):
        x = robjs.IntVector(x)
    else:
        x = robjs.FloatVector(x)

    ## now x is an array of numbers
    x = robjs.r.array(x, dim = shp, dimnames = dn)
    return x
    
def __enc2rds__(fi, fo):
    """ append encoded surface to R data """
    import rpy2
    import rpy2.robjects as robjs
    import rpy2.rlike.container as rlc
    R = rpy2.robjects.r
    from collections import OrderedDict

    tsk = hlp.load_pgz(fi)
    
    rbj = robjs.ListVector(OrderedDict())
    sbj = robjs.StrVector(tsk['sb'])
    vtx = robjs.StrVector(['V{:04X}'.format(i) for i in tsk['vi']])
    rbj.rx2['sbj'] = sbj
    rbj.rx2['vtx'] = vtx
    rbj.rx2['wms'] = robjs.StrVector([tsk['nm']])

    ## append surface encoding
    ## 2.enc: the encoding
    assert(tsk.has_key('enc'))
    
    od = OrderedDict()
    for k, v in tsk['enc'].iteritems():
        dn = robjs.ListVector(OrderedDict())
        dn.rx2['sbj']=robjs.StrVector(tsk['sb'])
        dn.rx2['vcd']=robjs.StrVector(
            ['C{:04X}'.format(i) for i in xrange(v.shape[1])])
        od['{}.{}'.format(*k)]=__np2r__(v, dn)

    rbj.rx2['enc'] = robjs.ListVector(od)
    
    ## append encoding to R data and save
    R.saveRDS(rbj, fo)
    del rbj, tsk, od
    print 'saved:', fo; sys.stdout.flush()
    print "xt: success"
    
def main(src, ovr = 0):
    from glob import glob as gg
    for fi in gg(pt.join(src, '*.pgz')):
        fo = pt.splitext(fi)[0] + '.rds'
        if pt.exists(fo) and not ovr:
            print 'exists:', fo
        else:
            __enc2rds__(fi, fo)
    
if __name__ == '__main__':
    import os
    import sys
    if len(sys.argv) > 1:
        wms = sys.argv[1]
    else:
        wms = None
        
    if len(sys.argv) > 2:
        src = sys.argv[2]
    else:
        src = pt.expandvars('$AZ_AENC')

    if wms is not None:
        main(src)
