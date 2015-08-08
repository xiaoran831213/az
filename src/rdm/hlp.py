import os
import numpy as np
import theano
import theano.tensor as T

def FX(fx = None):
    if fx is None:
        return theano.config.floatX
    else:
        theano.config.floatX = fx

## by default use 32bit float
FX('float32')

from theano import function as F

## * -------- random number helpers -------- * ##
rs_np = None                          # numpy random stream
rs_tn = None                          # theano random stream
__seed__ = 120
def set_seed(seed):
    global rs_np, rs_tn, __seed__
    rs_np = np.random.RandomState(seed)
    rs_tn = theano.tensor.shared_randomstreams.RandomStreams(
        rs_np.randint(2 ** 30))
    __seed__ = seed
    
def S(v, name = None, strict = False):
    """ create shared variable from v """
    ## wrap python type to numpy type
    if not isinstance(v, np.ndarray):
        v = np.array(v)

    ## wrap float type to default theano configuration
    if v.dtype is np.dtype('f8') and FX() is 'float32':
        v = np.asarray(v, dtype = 'f4')

    if v.dtype is np.dtype('i8') and FX() is 'float32':
        v = np.asarray(v, dtype = 'i4')

    if v.dtype is np.dtype('u8') and FX() is 'float32':
        v = np.asarray(v, dtype = 'u4')

    return theano.shared(v, name = name, strict = strict)

## rescaled values to [0, 1]
def rescale01(x, axis = None):
    """ rescale to [0, 1] """
    return (x - x.min(axis))/(x.max(axis) - x.min(axis))

## type checkers
def is_tvar(x):
    """ see if x is an theano symbolic is_tvar with
    no explict value. """
    return type(x) is T.TensorVariable

def is_tshr(x):
    """ show if a is_tvar is theano shared is_tvar """
    return type(x) is T.sharedvar.TensorSharedVariable

def is_tcns(x):
    """ show if a is_tvar is theano tensor constant """
    return type(x) is T.TensorConstant

def is_tnsr(x):
    """ see if x is an theano tensor """
    return (is_tvar(x) or is_tshr(x) or is_tcns(x))

## fetch parameters
def parms(y, chk = None):
    """
    find parameters in symbolic expression {y}. 

    chk: checker for allagible parameter. By default only shared
    variables could pass.
    """
    chk = is_tshr if chk is None else chk

    from collections import OrderedDict
    
    d = OrderedDict()
    q = [y]
    while len(q) > 0:
        v = q.pop()
        q.extend(v.get_parents())
        if chk(v):
            d[v] = v

    return d.keys()
            
# def save_pgz(fo, s):
#     """ save python object to gziped pickle """
#     import gzip
#     import cPickle
#     with gzip.open(fo, 'wb') as gz:
#         cPickle.dump(s, gz, cPickle.HIGHEST_PROTOCOL)

# def load_pgz(fi):
#     """ load python object from gziped pickle """
#     import gzip
#     import cPickle
#     with gzip.open(fi, 'rb') as gz:
#         return cPickle.load(gz)
    
