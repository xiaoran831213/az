import os.path as pt
import numpy as np
import sys
import cPickle

if not sys.path.count(pt.abspath('..')):
    sys.path.insert(0, pt.abspath('..'))
import hlp

def make_data():
    x0 = hlp.get_pk('dat/d48/1003', 0)
    x1 = hlp.get_pk('dat/d48/1035', 1)
    x2 = hlp.get_pk('dat/d48/2003', 2)
    x3 = hlp.get_pk('dat/d48/2035', 3)
    
    x = np.vstack((x0, x1, x2, x3))
    
    x = x.reshape(x.shape[0], -1)
    with open('dat/d48/full', 'wb') as pk:
        cPickle.dump(x, pk, cPickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    reload(hlp)
    pass
