import os
import numpy as np
import theano
import theano.tensor as T

def wrap_random():
    from theano.tensor.shared_randomstreams import RandomStreams as RS
    pass

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
            
        ## try evaluate v if it is a symbolic variable
        if is_symbol(v):
            try:
                val = v.eval()
            except theano.gof.fg.MissingInputError:
                val = None
            v = shared(val, name = v.name)
            ret.append(v)
            continue
        
        ## wrap python type to numpy type
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