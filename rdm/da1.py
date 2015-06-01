import os
import sys
import time

import numpy as np

import theano
import theano.tensor as T
from theano.tensor.shared_randomstreams import RandomStreams

from utils import tile_raster_images
import loader
import hlp

try:
    import PIL.Image as Image
except ImportError:
    import Image

FT = theano.config.floatX

class DA(object):
    """Denoising Auto-Encoder class (DA)

    """

    def __init__(
        self,
        np_rng,
        th_rng = None,
        n_vis =784,
        n_hid =500,
        t_w = None,
        t_b_hid = None,
        t_b_vis = None
    ):
        """
        Initialize the DA class by specifying the number of visible units (the
        dimension d of the input ), the number of hidden units ( the dimension
        d' of the latent or hidden space ) and the corruption level. The
        constructor also receives symbolic variables for the input, weights and
        bias. Such a symbolic variables are useful when, for example the input
        is the result of some computations, or when weights are shared between
        the DA and an MLP layer. When dealing with SdAs this always happens,
        the DA on layer 2 gets as input the output of the DA on layer 1,
        and the weights of the DA are used in the second stage of training
        to construct an MLP.

        """
        self.n_v = n_vis
        self.n_h = n_hid

        # create a Theano random generator that gives symbolic random values
        if not th_rng:
            th_rng = RandomStreams(np_rng.randint(2 ** 30))

        # note : W' was written as `W_prime` and b' as `b_prime`
        if not t_w:
            # W is initialized with `initial_W` which is uniformely sampled
            # from -4*sqrt(6./(n_vis+n_hid)) and
            # 4*sqrt(6./(n_hid+n_vis))the output of uniform if
            # converted using asarray to dtype
            # theano.config.floatX so that the code is runable on GPU
            initial_W = np.asarray(
                np_rng.uniform(
                    low=-4 * np.sqrt(6. / (n_hid + n_vis)),
                    high=4 * np.sqrt(6. / (n_hid + n_vis)),
                    size=(n_vis, n_hid)),
                dtype=theano.config.floatX)
            t_w = theano.shared(value=initial_W, name='W', borrow=True)

        if not t_b_vis:
            t_b_vis = theano.shared(value = np.zeros(n_vis, dtype = FT),
                name = 'b\'', borrow = True)

        if not t_b_hid:
            t_b_hid = theano.shared(value = np.zeros(n_hid, dtype = FT),
                name='b', borrow=True)

        self.t_w = t_w
        # b corresponds to the bias of the hidden
        self.t_b = t_b_hid
        # b_prime corresponds to the bias of the visible
        self.t_b_prime = t_b_vis
        # tied weights, therefore W_prime is W transpose
        self.t_w_prime = self.t_w.T
        self.th_rng = th_rng

        self.T_parm = [self.t_w, self.t_b, self.t_b_prime]

    def t_corrupt(self, T_X, T_lv):
        """This function keeps ``1-corruption_level`` entries of the inputs the
        same and zero-out randomly selected subset of size ``coruption_level``
        Note : first argument of theano.rng.binomial is the shape(size) of
               random numbers that it should produce
               second argument is the number of trials
               third argument is the probability of success of any trial

                this will produce an array of 0s and 1s where 1 has a
                probability of 1 - ``corruption_level`` and 0 with
                ``corruption_level``

                The binomial function return int64 data type by
                default.  int64 multiplicated by the T_X
                type(floatX) always return float64.  To keep all data
                in floatX when floatX is float32, we set the dtype of
                the binomial to floatX. As in our case the value of
                the binomial is always 0 or 1, this don't change the
                result. This is needed to allow the gpu to work
                correctly as it only support float32 for now.

        """
        return self.th_rng.binomial(
            size = T_X.shape, n = 1,
            p = 1 - T_lv,
            dtype = FT) * T_X

    def t_encode(self, T_X):
        """ Computes the values of the hidden layer """
        return T.nnet.sigmoid(T.dot(T_X, self.t_w) + self.t_b)

    def t_decode(self, T_X):
        """Computes the reconstructed input given the values of the
        hidden layer

        """
        return T.nnet.sigmoid(T.dot(T_X, self.t_w_prime) + self.t_b_prime)

    def f_train(self, t_data, t_corrupt = 0.2, t_rate = 0.1):
        """ return training function of the following signiture:
        input:
            lower and upper indices on training data
            alternative training data
        return:
            likelihood based cost
            square distance between training data and prediction
        
        """

        x = T.matrix('x')     # pipe data through this symble
        q = self.t_corrupt(x, t_corrupt)
        h = self.t_encode(q)
        z = self.t_decode(h)

        L = - T.sum(x * T.log(z) + (1 - x) * T.log(1 - z), axis=1)
        cost = T.mean(L)    # to be returned

        dist = T.mean(T.sqrt(T.sum((x - z) ** 2, axis = 1)))    # to be returned

        grad = T.grad(cost, self.T_parm)

        diff = [(p, p - t_rate * g) for p, g in zip(self.T_parm, grad)]

        t_fr = T.iscalar()
        t_to = T.iscalar()
        t_batch
        return theano.function(
            [t_batch],
            [cost, dist],
            updates = diff,
            givens = {x : t_data[t_batch]},
            name = "DA_trainer")
        
    def F_predictor(self):
        T_x = T.matrix('x')
        T_h = self.t_encode(T_x)
        T_z = self.t_decode(T_h)
        F_pred = theano.function(
            [T_x],
            T_z)
        return F_pred

def test_2(learning_rate = 0.1, training_epochs = 15,
            batch_size=20, output_folder='dA_plots'):

    import cPickle
    with open('dat/d48/2035') as pk:
        x = cPickle.load(pk)
        x = x.reshape(x.shape[0], -1)
        y = np.full(x.shape[0], 0)
    
    x = np.asarray(x, dtype = FT)
    y = np.asarray(y, dtype = FT)

    S_x = theano.shared(x, borrow = True)
    S_y = theano.shared(y, borrow = True)
    
    # compute number of minibatches for training, validation and testing
    s_batch = batch_size
    n_batch = S_x.get_value(borrow=True).shape[0] / s_batch

    hlp.mk_dir(output_folder)

    #####################################
    # BUILDING THE MODEL CORRUPTION 30% #
    #####################################
    np_rng = np.random.RandomState(123)
    th_rng = RandomStreams(np_rng.randint(2 ** 30))

    da = DA(
        np_rng = np_rng,
        th_rng = th_rng,
        n_vis = 48**3,
        n_hid = 100)

    train = da.f_train(t_data = S_x, t_corrupt = 0.2, t_rate = 0.1)

    ## -------- TRAINING --------
    start_time = time.clock()
    # go through training epochs
    for epoch in xrange(training_epochs):
        # go through trainng set
        c, d = [], []                     # cost, dist
        for i_batch in xrange(n_batch):
            r = train(i_batch * s_batch, (i_batch + 1) * s_batch)
            c.append(r[0])
            d.append(r[1])
        print 'Training epoch %d, cost %f, dist %f' % (epoch, np.mean(c), np.mean(d))

    end_time = time.clock()
    training_time = (end_time - start_time)
    print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))

    # start-snippet-4
    image = Image.fromarray(tile_raster_images(
        X=da.t_w.get_value(borrow=True).T,
        img_shape=(48*6, 48*8), tile_shape=(12, 16),
        tile_spacing=(1, 1)))
    image.save('filters_corruption_30.png')
    # end-snippet-4
    return da

def auc_da(da, x):
    from sklearn.metrics import roc_auc_score
    x = x.reshape(x.shape[0], -1)
    z = da.F_predictor()(x)
    s = np.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s.mean()

def make_da():
    np_rng = np.random.RandomState(123)
    th_rng = RandomStreams(np_rng.randint(2 ** 30))

    da = DA(
        np_rng = np_rng,
        th_rng = th_rng,
        n_vis = 48**3,
        n_hid = 100)
    return da

if __name__ == '__main__':
    pass
