import os
import numpy as np
import theano
import theano.tensor as T

def wrap_random():
    from theano.tensor.shared_randomstreams import RandomStreams as RS
    pass

def wrap_shared(*variables):
    """
    Try to wrap variables into theano shared variable.
    Return the intact variable if it is a symbolic tensor
    or is already a shared variable
    """
    from theano import shared
    ret = []
    for v in variables:
        ## do nothing to NoneType
        if v is None:
            ret.append(v)
            continue;
            
        ## do nothing if v is a symbolic tensor or
        ## is already a shared variable
        if type(v).__name__.startswith('Tensor'):
            ret.append(v)
            continue
        
        ## first wrap python type to numpy type
        if not isinstance(v, np.ndarray):
            v = np.array(v)
        dt = v.dtype

        ## make sure to use theano float type
        if dt.name.startswith('float'):
            dt = theano.config.floatX

        v = shared(np.asarray(v, dtype = dt), borrow = True)
        ret.append(v)

    if len(ret) == 1:
        ret = ret[0]
    return ret

def cross_entrophy(x, y, axis = None):
    """ symbolic expression of cross entrophy """
    x, y = wrap_shared(x, y)
    total = -T.sum(x * T.log(y) + (1 - x) * T.log(1 - y), axis = axis) 
    return T.mean(total)

def square_l2_norm(x, y, axis = None):
    """ symbolic expression of squared L2 norm """
    x, y = wrap_shared(x, y)
    total = T.sum((x - y) ** 2, axis = axis)
    return T.mean(total)

def rescale01(x, axis = None):
    """ rescale to [0, 1] """
    return (x - x.min(axis))/(x.max(axis) - x.min(axis))










