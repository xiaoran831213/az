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

class dA(object):
    """Denoising Auto-Encoder class (dA)

    A denoising autoencoders tries to reconstruct the input from a corrupted
    version of it by projecting it first in a latent space and reprojecting
    it afterwards back in the input space. Please refer to Vincent et al.,2008
    for more details. If T_x is the input then equation (1) computes a partially
    destroyed version of T_x by means of a stochastic mapping q_D. Equation (2)
    computes the projection of the input into the latent space. Equation (3)
    computes the reconstruction of the input, while equation (4) computes the
    reconstruction error.

    .. math::

        \tilde{T_x} ~ q_D(\tilde{T_x}|T_x)                                     (1)

        y = s(T_W \tilde{T_x} + T_b)                                           (2)

        T_x = s(T_W' y  + T_b')                                                (3)

        L(T_x,z) = -sum_{k=1}^d [x_k \log z_k + (1-x_k) \log( 1-z_k)]      (4)

    """

    def __init__(
        self,
        np_rng,
        th_rng = None,
        T_input = None,
        n_vis =784,
        n_hid =500,
        T_W = None,
        T_bhid = None,
        T_bvis = None
    ):
        """
        Initialize the dA class by specifying the number of visible units (the
        dimension d of the input ), the number of hidden units ( the dimension
        d' of the latent or hidden space ) and the corruption level. The
        constructor also receives symbolic variables for the input, weights and
        bias. Such a symbolic variables are useful when, for example the input
        is the result of some computations, or when weights are shared between
        the dA and an MLP layer. When dealing with SdAs this always happens,
        the dA on layer 2 gets as input the output of the dA on layer 1,
        and the weights of the dA are used in the second stage of training
        to construct an MLP.

        """
        self.n_vis = n_vis
        self.n_hid = n_hid

        # create a Theano random generator that gives symbolic random values
        if not th_rng:
            th_rng = RandomStreams(np_rng.randint(2 ** 30))

        # note : W' was written as `W_prime` and b' as `b_prime`
        if not T_W:
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
            T_W = theano.shared(value=initial_W, name='W', borrow=True)

        if not T_bvis:
            T_bvis = theano.shared(value = np.zeros(n_vis, dtype = FT),
                name = 'b\'', borrow = True)

        if not T_bhid:
            T_bhid = theano.shared(value = np.zeros(n_hid, dtype = FT),
                name='b', borrow=True)

        self.T_W = T_W
        # b corresponds to the bias of the hidden
        self.T_b = T_bhid
        # b_prime corresponds to the bias of the visible
        self.T_b_prime = T_bvis
        # tied weights, therefore W_prime is W transpose
        self.T_W_prime = self.T_W.T
        self.theano_rng = th_rng
        # if no input is given, generate a variable representing the input
        if T_input is None:
            # we use a matrix because we expect a minibatch of several
            # examples, each example being a row
            self.T_x = T.dmatrix(name='input')
        else:
            self.T_x = T_input

        self.T_parm = [self.T_W, self.T_b, self.T_b_prime]

    def get_corrupted_input(self, input, corruption_level):
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
                default.  int64 multiplicated by the input
                type(floatX) always return float64.  To keep all data
                in floatX when floatX is float32, we set the dtype of
                the binomial to floatX. As in our case the value of
                the binomial is always 0 or 1, this don't change the
                result. This is needed to allow the gpu to work
                correctly as it only support float32 for now.

        """
        return self.theano_rng.binomial(
            size = input.shape, n=1,
            p=1 - corruption_level,
            dtype = FT) * input

    def get_encoder(self, input):
        """ Computes the values of the hidden layer """
        return T.nnet.sigmoid(T.dot(input, self.T_W) + self.T_b)

    def get_decoder(self, hidden):
        """Computes the reconstructed input given the values of the
        hidden layer

        """
        return T.nnet.sigmoid(T.dot(hidden, self.T_W_prime) + self.T_b_prime)

    def get_cost_updates(self, corruption_level, learning_rate):
        """ This function computes the cost and the updates for one trainng
        step of the dA """

        tilde_x = self.get_corrupted_input(self.T_x, corruption_level)
        y = self.get_encoder(tilde_x)
        z = self.get_decoder(y)
        
        # note : we sum over the size of a datapoint; if we are using
        #        minibatches, L will be a vector, with one entry per
        #        example in minibatch
        T_L = - T.sum(self.T_x * T.log(z) + (1 - self.T_x) * T.log(1 - z), axis=1)
        
        # note : L is now a vector, where each element is the
        #        cross-entropy cost of the reconstruction of the
        #        corresponding example of the minibatch. We need to
        #        compute the average of all these to get the cost of
        #        the minibatch
        T_cost = T.mean(T_L)

        T_D = T.sqrt(T.sum((self.T_x - z) ** 2, axis = 1))
        T_dist = T.mean(T_D)

        # compute the gradients of the cost of the `dA` with respect
        # to its parameters
        T_grad = T.grad(T_cost, self.T_parm)

        # generate the list of updates
        T_updates = [
            (T_p, T_p - learning_rate * T_g)
            for T_p, T_g in zip(self.T_parm, T_grad)]

        return (T_cost, T_dist, T_updates)

    def get_predictor(self):
        T_x = T.matrix('x')
        T_h = self.get_encoder(T_x)
        T_z = self.get_decoder(T_h)
        F_pred = theano.function(
            [T_x],
            T_z)
        return F_pred
    
def test_dA(learning_rate = 0.1, training_epochs = 15,
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

    # start-snippet-2
    # allocate symbolic variables for the data
    T_i = T.lscalar('index')    # T_i to a [mini]batch
    T_x = T.matrix('x')  # the data is presented as rasterized images

    hlp.mk_dir(output_folder)

    #####################################
    # BUILDING THE MODEL CORRUPTION 30% #
    #####################################
    np_rng = np.random.RandomState(123)
    th_rng = RandomStreams(np_rng.randint(2 ** 30))

    da = dA(
        np_rng = np_rng,
        th_rng = th_rng,
        T_input = T_x,
        n_vis = 48**3,
        n_hid = 100)

    T_cost, T_dist, T_updates = da.get_cost_updates(
        corruption_level=0.2,
        learning_rate=learning_rate)

    train_da = theano.function(
        [T_i],
        [T_cost, T_dist],
        updates = T_updates,
        givens = {
            T_x: S_x[T_i * batch_size: (T_i + 1) * batch_size]
        })

    ## -------- TRAINING --------
    start_time = time.clock()
    # go through training epochs
    for epoch in xrange(training_epochs):
        # go through trainng set
        c, d = [], []                     # cost, dist
        for i_batch in xrange(n_batch):
            r = train_da(i_batch)
            c.append(r[0])
            d.append(r[1])
        print 'Training epoch %d, cost %f, dist %f' % (epoch, np.mean(c), np.mean(d))
    end_time = time.clock()

    training_time = (end_time - start_time)

    print >> sys.stderr, ('ran for %.2fm' % (training_time / 60.))
    # end-snippet-3

    # start-snippet-4
    image = Image.fromarray(tile_raster_images(
        X=da.T_W.get_value(borrow=True).T,
        img_shape=(48*6, 48*8), tile_shape=(12, 16),
        tile_spacing=(1, 1)))
    image.save('filters_corruption_30.png')
    # end-snippet-4

    return da

def auc_da(da, x):
    from sklearn.metrics import roc_auc_score
    z = da.get_predictor()(x)
    s = np.array([roc_auc_score(x[i], z[i]) for i in xrange(x.shape[0])])
    return s
    
if __name__ == '__main__':
    da = test_dA()
