import numpy as np
import nnt
from nnt import Nnt
import hlp
from hlp import T
from hlp import S
import pdb

class SA1(Nnt):
    """
    Stacked autoencoder
    """
    def __init__(self, dim):
        """
        Initialize the stacked auto encoder by a list of code dimensions.
        The weight and bias terms in the AEs are initialized by default rule.
        
        -------- parameters --------
        dim: a list of code dimensions
        d_1 is the dimension of the input (visible units)
        d_N is the dimension of the final code
        """
        super(SA1, self).__init__()
        
        import ae
        import cat
        
        """ the default view of the stack """
        sa = [ae.AE(d) for d in zip(dim[:-1], dim[1:])]
        
        """ the encoder view of the stacked autoencoder """
        ec = cat.Cat([a.ec for a in sa])

        """ the decoder view of the stacked autoencoder """
        dc = cat.Cat([a.dc for a in reversed(sa)])

        self.sa = sa                      # default view
        self.extend(sa)
        ## dimension of the stack decoder
        self.dim = dim
        self.dim.extend(reversed(dim[:-1]))

        self.ec = ec                      # encoder view
        self.dc = dc                      # decoder view

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

def test_sa1():
    import os.path as pt
    hlp.set_seed(120)

    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2, d/4]
    m = SA1(dim=dim)
    return x, m

if __name__ == '__main__':
    pass
