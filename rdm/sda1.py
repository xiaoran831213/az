import os
import sys
import time

import numpy

import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams

#from logistic_sgd import LogisticRegression, load_data
import da1

# start-snippet-1
class SDA(object):
    """Stacked denoising auto-encoder class (SdA)

    """

    def __init__(
        self,
        n_vis,
        np_rng,
        th_rng=None):

        """ This class is made to support a variable number of layers.

        """

        self.l_da = []
        self.parm = []

        if not th_rng:
            th_rng = RandomStreams(np_rng.randint(2 ** 30))

        self.n_vis = n_vis                # number of visible feature
        self.np_rng = np_rng
        self.th_rng = th_rng
        

    def add_DA(self, n_hid):
        if len(self.l_da):
            n_vis = self.l_da[-1].n_vis
        else:
            n_vis = self.n_vis

        da = da1.DA(
            np_rng = self.np_rng,
            th_rng = self.th_rng,
            n_vis = n_vis,
            n_hid = n_hid)
        
        self.l_da.append(da)
        self.parm.append(da.parm)
        
    def t_encode(self, t_x, layer = None):
        if not layer:
            layer = len(self.l_da)
        t_h = t_x
        for da in self.l_da[:layer]:
            t_h = da.t_encode(t_h)
        return t_h

    def t_decode(self, t_x, layer = None):
        pass

def test_sda():
    np_rng = np.random.RandomState(123)
    th_rng = RandomStreams(np_rng.randint(2 ** 30))

    sda = SDA(n_vis = 48**3, np_rng = np_rng)
    sda.add_DA(100)
    sda.add_DA(50)
    sda.add_DA(10)

    t_x = T.matrix('x')
    t_h = sda.t_encode(t_x)

def test_op(t_x):
    t_y = t_x + 1
    return t_y
    
def test_stack(t_x):
    h = t_x
    for i in xrange(3):
        h = test_op(h)
    return h

if __name__ == '__main__':
    test_sda()
