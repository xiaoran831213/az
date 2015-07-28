import os
import sys
import time
import os.path as pt
import numpy as np

import theano
import theano.tensor as T
from theano import shared as S
from theano.tensor.shared_randomstreams import RandomStreams
import hlp
import pdb

def __unref__(v):
    while callable(v):
        v = v()
    return v

def __isref__(v):
    return callable(v)

def __trans__(v, t = None):
    if t is None:
        return v
    if __isref__(v):
        return lambda :t(__unref__(v))
    if hlp.is_shared(v):
        return S(t(v).eval())
    return t(v)

class Lyr(object):
    """
    Denoising Auto-Encoder class (Lyr)
    """
    def __init__(
        self,
        d,
        np_rnd = None,
        th_rnd = None,
        w = None,
        b = None,
        s = None,
        x = None,
        tag = None
    ):
        """
        Initialize the neural network layer class by specifying the the dimension of the
        input, and the dimension of the output.
        The constructor also receives symbolic variables for the input, weights and bias.
        Such a symbolic variables are useful when, for example the input is the result of
        some computations, or when weights are shared between the layers
        
        -------- parameters --------
        d: a 2-tuple of input/output dimensions

        np_rnd: seed for Numpy random stream
        th_rnd: seed for Theano random stream

        w: (optional) weight of dimension (d_1, d_2), which is randomly filled by default.
        d_1 specify the input dimension
        d_2 specify the output dimension

        b: (optional) bias of dimension d_2, it is zero filled by default

        s: (optional) nonlinear tranformation of the weighted sum.
        By default the sigmoid function is adopted.
        To suppress nonlinearity, specify 1 instead.

        x: (optional) input data of dimension (N,d_1), here N is sample size.
        d_1 is the input dimension of the Lyr.
        """
        FT = theano.config.floatX

        ## I/O dimensions
        self.d = d

        if np_rnd is None:
            np_rnd = np.random.RandomState(120)

        ## create a Theano random generator that gives symbolic random values
        if th_rnd is None:
            th_rnd = RandomStreams(np_rnd.randint(2 ** 30))

        # note : W' was written as `W_prime` and b' as `b_prime`
        """
        # W is initialized with `initial_W` which is uniformely sampled
        # from -4*sqrt(6./(n_vis+n_hid)) and
        # 4*sqrt(6./(n_hid+n_vis))the output of uniform if
        # converted using asarray to dtype
        # theano.config.floatX so that the code is runable on GPU
        """
        if not w:
            initial_W = np.asarray(
                np_rnd.uniform(
                    low=-4 * np.sqrt(6. / (d[0] + d[1])),
                    high=4 * np.sqrt(6. / (d[0] + d[1])),
                    size=d),
                dtype=FT)
            w = theano.shared(value=initial_W, name='w', borrow=True)

        if not b:
            b = theano.shared(value = np.zeros(d[1], dtype = FT),
                name='b', borrow=True)

        if s is None:
            s = T.nnet.sigmoid
        self.s = s

        self.th_rnd = th_rnd
        
        self.__w__ = None       # weight term
        self.__v__ = None       # transformation on weight
        if w is not None:       # either assign or link the weight
            self.w(w)

        self.__b__ = None
        if b is not None:
            self.b(b)

        self.__x__ = T.matrix('x')
        if x is not None:
            self.x(x)

        self.tag = "" if tag is None else tag

    ## a Lyr cab be represented by the nonlinear funciton and I/O dimensions
    def __repr__(self):
        return '{}{}({}-{})'.format(
            self.tag, str(self.s)[0], self.d[0], self.d[1])

    def x(self, x = None):
        """
        getter/setter of layer input {x}
        -------- parameters --------
        x: the input object. could be one of the following:
        1) the output {y} from a lower layer, so as to wired the two layers;
        2)a Numpy array or compatible object, which turns the layer into an
        entry point of raw input.
        
        when {x} is None, the function serves as a getter of layer input
        returns:
        Theano Tensor variable if the layer is an intermidiate of a network
        Theano shared variable if the layer is an entry point of data
        """
        ## receive input through this symble
        if x is None:
            return self.__x__() if callable(self.__x__) else self.__x__
        else:
            if callable(x):
                self.__x__ = x
            else:
                if hlp.is_tensor(x):
                    self.__x__ = x
                else:
                    self.__x__ = S(x, 'x')

    def y(self):
        """
        getter of layer output {y}
        returns:
        Theano Tensor variable of y
        """
        x = self.x()
        w = self.w()
        b = self.b()
        y = T.dot(x, w) + b
        if self.s is 1:
            return  T.dot(x, w) + b
        else:
            return  self.s(T.dot(x, w) + b)

    def w(self, w = None, t = None):
        """
        getter/setter of weight {w}
        w: use shared tensor to assign independent parameter;
        use callable object to refer assign dependent parameter
        
        t: transformation over dependent parameter
        """
        if w is None:
            return __unref__(self.__w__)
        else:
            self.__w__ = __trans__(w, t)

    def b(self, b = None, t = None):
        """
        getter/setter of bias {b}
        """
        if b is None:
            return __unref__(self.__b__)
        else:
            self.__b__ = __trans__(b, t)

    def p(self):
        """
        getter of parameters {p}
        only list independent parameters
        """
        return [p for p in [w(), b()] if hlp.is_shared(p)]

    def itr_back(self, stop = None):
        """ the backward layer iterator, starts with the caller.
        
        stop: the stoping layer, the one below the last layer to visit.
        If unspecified, all layers will be through
        """
        l = self                # starts with the caller
        ## callable(l.__x__) being True means a lower layer is
        ## connected to this layer
        while l is not stop:
            r = l                 # remember current layer

            ## point to lower layer
            if callable(l.__x__):
                l = l.__x__.__self__
            else:
                stop = l
                
            ## yield current layer
            yield r

    def iter_p(self):
        """ iterate downwards to collect parameters """
        for l in self.itr_back():
            for p in l.p():
                yield p
            
def test_lyr():
    x = np.load(pt.expandvars('$AZ_IMG1/lh001F1.npz'))['vtx']['tck']
    ##x = np.asarray(x, dtype = T.config.floatX)
    x = x.reshape(x.shape[0], -1)
    d = (x.shape[1], x.shape[1]/2)
    x = (x - x.min()) / (x.max() - x.min())
    x = hlp.to_shared(x)
    
    np_rng = np.random.RandomState(120)
    nt = Lyr(d=d, np_rnd = np_rng, x=x)
    return nt

if __name__ == '__main__':
    theano.config.floatX = 'float32'
    pass
