import numpy as np
import nnt
from nnt import Nnt
import hlp
from hlp import T
from hlp import S
import pdb
import cat
from cat import Cat

class SAE(Cat):
    """
    Stacked Auto Encoder
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
        
        import ae
        
        """ the default view of the stacked autoencoders"""
        sa = [ae.AE(d) for d in zip(dim[:-1], dim[1:])]
        
        """ the encoder view of the stacked autoencoders """
        ec = cat.Cat([a.ec for a in sa])

        """ the decoder view of the stacked autoencoders """
        dc = cat.Cat([a.dc for a in reversed(sa)])

        self.sa = sa            # default view
        self.ec = ec            # encoder view
        self.dc = dc            # decoder view

        nnts = []
        nnts.extend(ec)
        nnts.extend(dc)
        super(SAE, self).__init__()
        
        ## dimension of the stack decoder
        # self.dim = dim
        # self.dim.extend(reversed(dim[:-1]))
        
def test_sa1():
    import os.path as pt
    hlp.set_seed(120)

    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2, d/4, d/8]
    m = SAE(dim=dim)
    return x, m

if __name__ == '__main__':
    pass










