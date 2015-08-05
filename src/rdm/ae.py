import numpy as np
from nnt import Nnt
import lyr
from lyr import Lyr
import hlp
from hlp import T
from hlp import S
import pdb

class AE(Nnt):
    """
    Basic auto-encoder
    """
    def __init__(
        self, dim,
        ec_w = None, ec_b = None, ec_s = None,
        dc_w = None, dc_b = None, dc_s = None):
        """
        Initialize the denosing auto encoder by specifying the the dimension of the input
        and output.
        The constructor also receives symbolic variables for the weights and bias.
        
        -------- parameters --------
        dim: a 2-tuple of input/output dimensions
        d_1 is the dimension of the input (visible units)
        d_2 is the dimension of the code (hidden units)

        ec_w: (optional) encoder weight of dimension (d_1, d_2),
        by default it is randomly initialized.

        ec_b: (optional) encoder bias of dimension d_2,
        by default the initial value is 0.

        ec_s: (optional) encoder nonlinearity. The default is Sigmoid.
        Specify 1 to disable nonlinearity
        
        dc_w: (optional) decoder weight of dimension (d_2, d_1), by default it is
        constrained to be the transpose of encoding weight

        dc_b: (optional) decoder bias of dimension d_1, by default the initial value is 0.
        dc_s: (optional) decoder nonlinearity, by default be the same with the encoder.
        """
        super(AE, self).__init__()
        self.dim = [dim[0], dim[1], dim[0]]

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
    import numpy as np
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = x.shape[1]
    x = hlp.rescale01(x)

    dim = [d/1, d/2]
    m = AE(dim=dim)

    return x, m

if __name__ == '__main__':
    pass
