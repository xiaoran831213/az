import numpy as np
import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams
from theano import pp
import hlp
from hlp import S
import pdb

def dist_ce(y, z):
    """ symbolic expression of cross entrophy
    y: produced output
    z: expected output

    The first dimension denote unit of sampling
    """
    d = - z * T.log(y) - (1 - z) * T.log(1 - y)
    bySample = T.sum(d.flatten(ndim = 2), axis = 1) 
    return T.mean(bySample)

def dist_l2(y, z):
    """ symbolic expression of L2 norm
    y: produced output
    z: expected output

    The first dimension denote unit of sampling
    """
    d = (y - z) ** 2
    bySample = T.sqrt(T.sum(d.flatten(ndim = 2), axis = 1))
    return T.mean(bySample)

class wreg_l2:
    """
    a callable class, returns symbolic expression of weight
    regulator in the form of L2 norm.
    """
    def __init__(self, lmb = None):
        """ the constructor.
        lmb: the relative importance of weight regulator.
        usually denoted as 'lambda' in literatures.
        """
        lmb = 0.1 if lmb is None else float(lmb)
        self.lmb = hlp.S(lmb)

    def __call__(self, lsW):
        """
        return the symbolic regulator term of a list of
        weights
        """
        w = T.concatenate([w.flatten() for w in lsW])
        return self.lmb * T.sqrt(T.sum(w ** 2))

    def __repr__(self):
        return 'f={}, l={}'.format('L2', self.lmb)

class Trainer(object):
    """
    Class for neural network training.
    """
    def __init__(
            self, nnt,
            src = None, xpt = None,
            call_dist=None,
            call_wreg=None,
            bsz=None, mmt=None, lrt=None):
        """
        Constructor.
        : -------- parameters -------- :
        nnt: the neural network to be trained

        src: training source, the first dimension of which stands for sample units.
        if unspecified, the trainer will try to evaluate the entry point and cache
        the result as source data.
        
        xpt: training expect, the first dimension should be identical with source.
        if unspecified, a self-supervied training is assumed and the expect is set
        to be identical to the source.

        call_dist: an expression builder for ouput-expect distance(or distortion).
        the supplied builder should return an expression given the symbolic output
        {y}, and the expect {z}. the expression must evaluate to a scalar.

        call_wreg: an expression builder for regulator term.
        The build should return an expression given the list of independant weights
        to be used to panalize. the expression must evaluate to a scalar
        for the regulator. restriction oevr the magnitude of weights supposingly
        helps the final modle more generalizable.

        bsz: size of training batch
        mmt: momentom
        lrt: learning rate
        lmb: weight decay speed (usually denoted with 'lambda' in literatures)
        """
        ## grand source and expect
        src = entry.eval() if src is None else src
        xpt = src if xpt is None else xpt
        self.src = S(src, 'src')
        self.xpt = S(xpt, 'xpt')
        
        ## form of loss function
        self.call_dist = dist_ce if call_dist is None else call_dist

        ## form of weight regulator
        self.call_wreg = (lambda w:S(0.0, 'lmb')) if call_wreg is None else call_wreg
        
        ## current epoch index
        self.eph = S(0, 'eph')
            
        ## training batch ppsize
        bsz = 20 if bsz is None else bsz
        self.bsz = S(bsz)

        ## current batch index
        self.bat = S(0)

        ## momentumn, make sure momentum is a sane value
        mmt = 0.0 if mmt is None else mmt
        assert mmt < 1 and mmt >= 0
        self.mmt = S(mmt, 'mmt')

        ## learning lrt
        lrt = 0.1 if lrt is None else lrt
        self.lrt = S(lrt, 'lrt')

        ## -------- construct trainer function -------- *
        ## 1) symbolic expressions
        x = T.matrix('x') # the symbolic batch source
        z = T.matrix('z') # the symbolic batch expect
        
        y = exitp()      # get symbolic batch output from network exit
        entry(old_wire)  # wire old source back to the network entry

        ## list of independant symbolic parameters to be tuned
        parm = list(exitp.__self__.itr_p())

        ## list of independant symbolic weights to apply decay
        lswt = [l.w() for l in exitp.__self__.itr_back() if hlp.is_tshr(l.w())]

        ## symbolic cost
        dist = self.call_dist(y, z)
        wreg = self.call_wreg(lswt)
        cost = dist + wreg

        ## symbolic gradient of cost WRT parameters
        grad = T.grad(cost, parm)
        gsum = T.sum([T.abs_(g).sum() for g in grad])

        ## 2) define updates after each batch training
        ZPG = zip(parm, grad)
        up = []
        ## update parameters using gradiant decent, and momentum
        for p, g in ZPG:
            ## initialize accumulated gradient
            h = S(np.zeros_like(p.eval()))

            ## accumulate gradient, partially historical (due to the momentum),
            ## partially noval
            up.append((h, self.mmt * h + (1 - self.mmt) * g))

            ## update parameters by stepping down the accumulated gradient
            up.append((p, p - self.lrt * h))

        ## update batch and eqoch index
        bnx = (((self.bat + 1) * self.bsz) % src.shape[0]) / self.bsz
        enx = self.eph + ((self.bat + 1) * self.bsz) / src.shape[0]
        up.append((self.bat, bnx))
        up.append((self.eph, enx))

        ## 3) the trainer functions
        from hlp import F
        ## feed symbols with explicit data in batches
        gvn = {
            x: self.src[self.bat * self.bsz : (self.bat + 1) * self.bsz],
            z: self.xpt[self.bat * self.bsz : (self.bat + 1) * self.bsz]}

        ## each invocation sends one batch of training examples to the network,
        ## calculate total cost and tuen the parameters using gradient decent.
        self.step = F([], cost, name = "step", givens = gvn, updates = up)

        ## return intermidiate result for batch training
        self.dist = F([], dist, name = "dist", givens = gvn)
        self.wreg = F([], wreg, name = "wreg")
        self.cost = F([], cost, name = "cost", givens = gvn)
        self.grad = dict([(p, F([], g, givens = gvn)) for p, g in ZPG])
        self.gsum = F([], gsum, name = "gsum", givens = gvn)
        ## -------- done with trainer functions -------- *

    def tune(self, nep = 1, npt = 1):
        """ tune the parameters by running the trainer {nep} epoch.
        an epoch is one going through of all samples
        """
        b0 = self.bat.get_value()
        e0 = self.eph.get_value()
        eN = e0 + nep
        pN = e0 + npt
        while self.eph.get_value() < eN or self.bat.get_value() < b0:
            self.step()
            i = self.eph.get_value().item()
            j = self.bat.get_value().item()
            d = self.dist().item()
            g = self.gsum().item()
            if  i < pN or j < b0:
                continue
            print 'e{:04d}: {:08.3f} {:07.4f}'.format(i, d, g)
            pN = i + npt

def data_x():
    import os.path as pt
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    x = x.reshape(x.shape[0], -1)
    d = (x.shape[1], x.shape[1]/2)
    x = (x - x.min()) / (x.max() - x.min())
    return x

def data_n(x):
    import lyr
    reload(lyr)
    from lyr import Lyr

    d = x.shape[1]
    rnd = np.random.RandomState(150)
    l1 = Lyr(dim = (d/1, d/2), np_rnd = rnd, tag = 'E1', x = x)
    l2 = Lyr(dim = (d/2, d/1), np_rnd = rnd, tag = 'D1', x = l1.y)
    l2.w(l1.w, T.transpose)
    return l1, l2

def test_material():
    x = data_x()
    n = data_n(x)
    t = Trainer(n[0].x, n[1].y, src = x, xpt = x, lrt = 0.01)
    return t

if __name__ == '__main__':
    theano.config.floatX = 'float32'
    pass
