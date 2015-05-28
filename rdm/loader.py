import os.path as pt
import numpy as np
import sys
import theano
import theano.tensor as T
import cPickle

if not sys.path.count(pt.abspath('..')):
    sys.path.insert(0, pt.abspath('..'))
import hlp

    
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

def make_data():
    x0 = hlp.get_pk('dat/d32', 0)
    y0 = np.full(x0.shape[0], 0)
    x1 = hlp.get_pk('dat/d32', 1)
    y1 = np.full(x1.shape[1], 1)
    x2 = hlp.get_pk('dat/d32', 2)
    y2 = np.full(x2.shape[2], 2)
    x3 = hlp.get_pk('dat/d32', 3)
    y3 = np.full(x3.shape[3], 3)
    
    x = np.vstack((x0, x1, x2, x3))
    y = np.hstack((y0, y1, y2, y3))
    
    x = x.reshape(x.shape[0], -1)
    with open('dat/tst/t32', 'wb') as pk:
        cPickle.dump((x, y), pk, cPickle.HIGHEST_PROTOCOL)
    return (x, y)

if __name__ == "__main__":
    pass
