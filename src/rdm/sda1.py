import numpy as np
import lyr
from lyr import Lyr
import sd1
from sd1 import DA1
import hlp
from hlp import T
from hlp import S
import pdb

class SA1(Nnt):
    """
    Generic layer of neural network
    """
    def __init__(self, dim):
        """
        Initialize the denosing auto encoder by specifying the the dimension of the input
        and output.
        The constructor also receives symbolic variables for the weights and bias.
        
        -------- parameters --------
        dim: a list of code dimensions
        d_1 is the dimension of the input (visible units)

        """
        super(SA1, self).__init__()
        self.dim = [dim[0], dim[1], dim[0]]

        ZDT = zip(dims[:-1], dims[1:], xrange(len(dims)))

        das = [DA1(d), for d in zip(dims[:-1], dims[1:])]
        """ the encoder view of the autoencoder """
        ec = Lyr(dim, ec_w, ec_b, ec_s)

        ## dimension of the decoder
        dim = [dim[1], dim[0]]
        
        ## ec_w prime, the decoding weight, remain None to enforce the dc_w = ec_w.T constraint
        dc_w = ec.w.T if dc_w is None else dc_w

        """ the decoder view of the autoencoder """
        dc_s = ec.s if dc_s is None else dc_s

        """ the decoder view of the DA """
        dc = Lyr(dim, dc_w, dc_b, dc_s)

        self.append(ec)
        self.append(dc)
        self.ec = ec
        self.dc = dc

    def y(self, x):
        """
        build symbolic expression of latent code given input {x}
        """
        return self.ec(x)

    def z(self, x):
        """
        build symbolic expression of reconstruction of input {x}
        """
        return self.dc(self.ec(x))

def test_da1():
    import os.path as pt

    hlp.set_seed(120)
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = (x.shape[1], x.shape[1]/2)
    x = hlp.rescale01(x)

    m = DA1(dim=d)
    return x, m

if __name__ == '__main__':
    pass
