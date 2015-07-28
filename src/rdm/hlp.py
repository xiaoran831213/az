import os
import numpy as np
import theano
import theano.tensor as T
FT = theano.config.floatX

def wrap_random():
    from theano.tensor.shared_randomstreams import RandomStreams as RS
    pass

def S(v, name = None, strict = False):
    """ create shared variable from v """

    ## wrap python type to numpy type
    if not isinstance(v, np.ndarray):
        v = np.array(v)

    ## wrap float type to default theano configuration
    if v.dtype in (np.float32, np.float64) and v.dtype is not FT:
        v = np.asarray(v, dtype = FT)

    return theano.shared(v, name = name, strict = strict)
    
def to_shared(*variables):
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

        ## do nothing if v is aready a shared variable
        if is_shared(v):
            ret.append(v)
            continue
            
        ## wrap python type to numpy type
        if not isinstance(v, np.ndarray):
            v = np.array(v)
        dt = v.dtype

        ## make sure to use theano float type
        if dt in (np.float32, np.float64) and dt is not FT:
            v = np.asarray(v, dtype = FT)
            
        v = shared(v, borrow = True)
        ret.append(v)

    if len(ret) == 1:
        ret = ret[0]
    return ret


def is_symbol(x):
    """ see if x is an theano symbolic variable with
    no explict value. """
    return type(x) is T.TensorVariable

def is_shared(x):
    """ show if a variable is theano shared variable """
    return type(x) is T.sharedvar.TensorSharedVariable

def is_tensor(x):
    """ see if x is an theano tensor """
    return (is_symbol(x) or is_shared(x))
    
def rescale01(x, axis = None):
    """ rescale to [0, 1] """
    return (x - x.min(axis))/(x.max(axis) - x.min(axis))
