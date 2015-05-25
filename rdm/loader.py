import os.path as pt
import numpy as np
import sys
import theano
import theano.tensor as T

if not sys.path.count(pt.abspath('..')):
    sys.path.insert(0, pt.abspath('..'))

import cPickle
import hlp

def mk_dat(src, fo, ovr = False):
    print "mk_dat: ", src, " -> ", fo
    renew = False
    if pt.isfile(fo):     # skip exists
        if not ovr:
            print fo, "exists"
            return
        else:
            renew = True

    x = []
    y = []
    for o in hlp.itr_pk(src):
        x.append(o['x'])
        y.append(o['y'])

    x = np.asarray(x, dtype = 'u1')
    y = np.asarray(y, dtype = 'u1')

    o = {'x': x, 'y': y}

    with open(fo, 'wb') as pk:
        cPickle.dump(o, pk, cPickle.HIGHEST_PROTOCOL)

    if renew:
        print fo, "renewed"
    else:
        print fo, "created"

def gt_dat(fi, fr = 0, to = 64):
    with open(fi, 'rb') as pk:
        d = cPickle.load(pk)
        
    y = d['y']
    x = d['x']
    print 'slice:', fr, to
    x = x[:, fr:to, fr:to, fr:to].reshape(
        x.shape[0], -1)

    shared_x = theano.shared(
        np.asarray(
            x, dtype = theano.config.floatX),
        borrow = True)

    shared_y = theano.shared(
        np.asarray(y, dtype = theano.config.floatX),
        borrow = True)
    
    return shared_x, shared_y

if __name__ == "__main__":
    pass
#    del hlp, pt, np
