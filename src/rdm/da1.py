import os.path as pt
import numpy as np
import lyr
from lyr import Lyr
import hlp
from hlp import T
from hlp import S
import pdb

class DA1(Lyr):
    """
    Generic layer of neural network
    """
    def __init__(
        self, dim, w = None, b = None, s = None,
        w_p = None, b_p = None, s_p = None, tag = None):
        """
        Initialize the denosing auto encoder by specifying the the dimension of the input
        and output.
        The constructor also receives symbolic variables for the weights and bias.
        
        -------- parameters --------
        dim: a 2-tuple of input/output dimensions

        w: (optional) weight of dimension (d_1, d_2), which is randomly filled by default.
        d_1 specify the input dimension
        d_2 specify the output dimension

        b: (optional) bias of dimension d_2, it is zero filled by default

        s: (optional) nonlinear tranformation of the weighted sum.
        By default the sigmoid function is used.
        To disable nonlinearity, specify 1 instead.
        
        w_p: (optional) reconstruction weight of dimension (d_2, d_1), by default it is
        constrained to be the transpose of encoding weight

        b_p: (optional) decoder bias of dimension d_1, filled with 0 by default
        s_p: (optional) decoder non-linear function, by default the same with the encoder
        """

        super(DA1, self).__init__(dim, w, b, s, tag)

        ## w prime, the decoding weight, remain None to enforce the w_p = w.T constraint
        self.w_p = self.w.T if w_p is None else w_p

        ## b prime, the decoding bias
        b_p = S(np.zeros(dim[0], dtype = hlp.FT), name='b_p') if b_p is None else b_p
        self.b_p = b_p

        ## s prime, the decoding non-linear transformation
        s_p = self.s if s_p is None else s_p
        self.s_p = s_p

        """ the encoder view of the DA """
        self.ec = Lyr(dim, self.w, self.b, self.s)
        """ the decoder view of the DA """
        self.dc = Lyr((dim[1], dim[0]), self.w_p, self.b_p, self.s_p)
        

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

    def p(self):
        """
        return independant parameters.
        """
        ret = []
        ret.extend(self.ec.p())
        ret.extend(self.dc.p())
        return ret

    def __repr__(self):
        return '{}{}-{}'.format(self.tag, self.ec, self.dc)


def test_da1():

    hlp.set_seed(120)
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    d = (x.shape[1], x.shape[1]/2)
    x = hlp.rescale01(x)

    m = DA1(dim=d)
    return x, m

if __name__ == '__main__':
    pass
