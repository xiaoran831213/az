import numpy as np
import theano
import hlp
from hlp import S
from hlp import T
import pdb
import nnt
from nnt import Nnt

class Serial(Nnt):
    """
    Generic serial structured neural network
    """
    def __init__(self, nnts):
        """
        Initialize the neural network layer class by specifying the the dimension of the
        input, and the dimension of the output.
        The constructor also receives symbolic variables for the input, weights and bias.
        Such a symbolic variables are useful when, for example the input is the result of
        some computations, or when weights are shared between the layers
        
        -------- parameters --------
        dim: a 2-tuple of input/output dimensions

        w: (optional) weight of dimension (d_1, d_2), which is randomly filled by default.
        d_1 specify the input dimension
        d_2 specify the output dimension

        b: (optional) bias of dimension d_2, it is zero filled by default

        s: (optional) nonlinear tranformation of the weighted sum.
        By default the sigmoid function is used.
        To suppress nonlinearity, specify 1 instead.
        """
        super(Serial, self).__init__()

        ## I/O dimensions
        self.dim = (nnts[0].dim[0], nnts[-1].dim[-1])

        self.extend(nnts)

    def y(self, x):
        """
        build sybolic expression of layer output {y} given input {x}
        this also the defaut expression returned when the Lyr object is
        called as a function
        """
        y = x
        for net in self:
            y = net.y(x)
        return y

def test_serial():
    from os import path as pt
    from lyr import Lyr
    hlp.set_seed(120)
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = (x.shape[1], x.shape[1]/2)
    x = hlp.rescale01(x)
    
    nt = Lyr(dim=d)
    return x, nt

if __name__ == '__main__':
    pass
