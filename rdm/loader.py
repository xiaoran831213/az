import os.path as pt
import numpy as np
import sys
import theano
import theano.tensor as T
import cPickle

if not sys.path.count(pt.abspath('..')):
    sys.path.insert(0, pt.abspath('..'))
import hlp

def make_data():
    x0 = hlp.get_pk('dat/d48', 0)
    y0 = np.full(x0.shape[0], 0)
    x1 = hlp.get_pk('dat/d48', 1)
    y1 = np.full(x1.shape[0], 1)
    x2 = hlp.get_pk('dat/d48', 2)
    y2 = np.full(x2.shape[0], 2)
    x3 = hlp.get_pk('dat/d48', 3)
    y3 = np.full(x3.shape[0], 3)
    
    x = np.vstack((x0, x1, x2, x3))
    y = np.hstack((y0, y1, y2, y3))
    
    x = x.reshape(x.shape[0], -1)
    with open('dat/tst/t48', 'wb') as pk:
        cPickle.dump((x, y), pk, cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    reload(hlp)
    pass
