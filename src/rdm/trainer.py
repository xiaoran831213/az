import os
import sys
import time
import os.path as pt
import numpy as np
import theano
import theano.tensor as T
from theano import function as F
from theano.tensor.shared_randomstreams import RandomStreams
from theano import pp
from theano import shared as S
import hlp
import pdb

def cross_entrophy(y, z):
    """ symbolic expression of cross entrophy
    y: produced output
    z: expected output

    The first dimension denote unit of sampling
    """
    d = - z * T.log(y) - (1 - z) * T.log(1 - y)
    total = T.sum(d.flatten(ndim = 2), axis = 1) 
    return T.mean(total)

def square_l2_norm(y, z):
    """ symbolic expression of squared L2 norm
    y: produced output
    z: expected output

    The first dimension denote unit of sampling
    """
    d = (y - z) ** 2
    total = T.sqrt(T.sum(d.flatten(ndim = 2), axis = 1))
    return T.mean(total)

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
        lmb = 0.1 if lmb is None else lmb
        self.lmb = hlp.to_shared(lmb)
        self.set_lambda = self.lmb.set_value
        self.get_lambda = self.lmb.get_value

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
            self,
            entry, exitp,
            src = None, xpt = None,
            call_dist=None,
            call_wreg=None,
            bsz=None, mmt=None, lrt=None):
        """
        Constructor.
        : -------- parameters -------- :
        entry: network entry point
        exitp: network exit
        entry and exitp must be connected

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
        self.call_dist = cross_entrophy if call_dist is None else call_dist

        ## form of weight regulator
        self.call_wreg = (lambda w:S(0.0, 'lmb')) if call_wreg is None else call_wreg
        
        ## current epoch index
        self.eph = S(0, 'eph')
            
        ## training batch ppsize
        bsz = 20 if bsz is None else bsz
        self.bsz = hlp.to_shared(bsz)

        ## current batch index
        self.bat = hlp.to_shared(0)

        ## momentumn, make sure momentum is a sane value
        mmt = 0 if mmt is None else mmt
        assert mmt < 1 and mmt >= 0
        self.mmt = S(mmt, 'mnt')

        ## learning lrt
        lrt = 0.1 if lrt is None else lrt
        self.lrt = S(lrt, 'lrt')

        ## -------- construct trainer function -------- *
        ## 1) symbolic expressions
        x = T.matrix('x') # the symbolic batch source
        z = T.matrix('z') # the symbolic batch expect
        
        old_wire = entry()      # backup
        entry(x)         # wire the batch source to the network entry
        y = exitp()      # get symbolic batch output from network exit
        entry(old_wire)  # wire old source back to the network entry

        ## list of independant symbolic parameters to be tuned
        parm = list(exitp.__self__.iter_p())

        ## list of independant symbolic weights to apply decay
        lswt = [l.w() for l in exitp.__self__.itr_back() if hlp.is_shared(l.w())]

        ## symbolic cost
        dist = self.call_dist(y, z)
        wreg = self.call_wreg(lswt)
        cost = dist + wreg

        ## symbolic gradient of cost WRT parameters
        grad = T.grad(cost, parm)
        pg = zip(parm, grad)

        ## 2) define updates after each batch training
        up = []
        ## update parameters using gradiant decent, and momentum
        for p, g in pg:
            ## initialize accumulated gradient history
            h = theano.shared(0 * p.get_value())

            ## update gradient accumulate, part history (the momentum),
            ## part new
            up.append((h, self.mmt *h + (1 - self.mmt) * g))

            ## update parameters by stepping down the accumulated gradient
            up.append((p, p - self.lrt * h))

        ## update batch and eqoch index
        bnx = (((self.bat + 1) * self.bsz) % src.shape[0]) / self.bsz
        enx = ((self.bat + 1) * self.bsz) / src.shape[0]
        up.append((self.bat, bnx))
        up.append((self.eph, enx))

        ## 3) the trainer functions
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
        self.grad = dict([(p, F([], g, givens = gvn)) for p, g in pg])
        ## -------- done with trainer functions -------- *

def data_x():
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
    rnd = np.random.RandomState(120)
    l1 = Lyr(d = (d/1, d/2), np_rnd = rnd, tag = 'E1', x = x)
    l2 = Lyr(d = (d/2, d/4), np_rnd = rnd, tag = 'E2', x = l1.y)
    l3 = Lyr(d = (d/4, d/8), np_rnd = rnd, tag = 'E3', x = l2.y)

    l4 = Lyr(d = (d/8, d/4), np_rnd = rnd, tag = 'D3', x = l3.y)
    l5 = Lyr(d = (d/4, d/2), np_rnd = rnd, tag = 'D2', x = l4.y)
    l6 = Lyr(d = (d/2, d/1), np_rnd = rnd, tag = 'D1', x = l5.y)

    return (l1, l2, l3, l4, l5, l6)

def test_material():
    x = data_x()
    nt1 = data_n(x)
    nt2 = data_n(x)
    t1 = Trainer(nt1[0].x, nt1[5].y, src = x, xpt = x, lrt = 0.1, mmt = 0.7)
    t2 = Trainer(nt2[0].x, nt2[5].y, src = x, xpt = x, lrt = 0.1)
    return(t1, t2)

def text_trainer():
    pass

if __name__ == '__main__':
    theano.config.floatX = 'float32'
    pass
