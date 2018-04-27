"""
Qiwen Hu 2018
Variation Autoencoder for single cell data
Modified from Greg Way's code https://github.com/greenelab/tybalt/blob/master/scripts/vae_pancancer.py#L14
scripts/vae_singlecell.py

Output:
1) Loss and validation loss for models
2) Models/Z-matrix for different parameters
3) t-sne on tybalt and/or rna-seq features

"""
import os
import argparse
import numpy as np
import pandas as pd

import tensorflow as tf
from keras.layers import Input, Dense, Lambda, Layer, Activation
from keras.layers.normalization import BatchNormalization
from keras.models import Model, Sequential
from keras import backend as K
from keras import metrics, optimizers
from keras.callbacks import Callback
from sklearn import manifold
import string

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--input_file',
                    help='input_file')
parser.add_argument('-d', '--depth',
                    help='Number of layers between input and latent layer')
parser.add_argument('-t', '--tsne_RNAseq', default=False,
                    help='Perform t-SNE on single cell gene expression features, True or False')
parser.add_argument('-n', '--num_components', default=20,
		    help='The latent space dimensionality to test')

args = parser.parse_args()

# Set hyper parameters
input_file = args.input_file
depth = int(args.depth)
tsne_RNAseq = args.tsne_RNAseq
latent_dim = int(args.num_components)

# Set parameters
first_layer = 100
second_layer = 100
batch_size = 50
epochs = 200
learning_rate = 0.0005
kappa = 1

# Load Data
rnaseq_file = os.path.join('data', input_file)
rnaseq_df = pd.read_table(rnaseq_file, index_col=0)
rnaseq_df = rnaseq_df.T

# Set architecture dimensions
original_dim = rnaseq_df.shape[1]
epsilon_std = 1.0
beta = K.variable(0)
if depth == 2:
    hidden_dim = int(first_layer)

if depth == 3:
    hidden_dim = int(first_layer)
    hidden_dim_2 = int(second_layer)

# Random seed
seed = int(np.random.randint(low=0, high=10000, size=1))
np.random.seed(seed)

# Function for reparameterization trick to make model differentiable
def sampling(args):

    # Function with args required for Keras Lambda function
    z_mean, z_log_var = args

    # Draw epsilon of the same shape from a standard normal distribution
    epsilon = K.random_normal(shape=tf.shape(z_mean), mean=0.,
                              stddev=epsilon_std)

    # The latent vector is non-deterministic and differentiable
    # in respect to z_mean and z_log_var
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z


class CustomVariationalLayer(Layer):
    """
    Define a custom layer that learns and performs the training

    """
    def __init__(self, **kwargs):
        # https://keras.io/layers/writing-your-own-keras-layers/
        self.is_placeholder = True
        super(CustomVariationalLayer, self).__init__(**kwargs)

    def vae_loss(self, x_input, x_decoded):
        reconstruction_loss = original_dim * \
                              metrics.binary_crossentropy(x_input, x_decoded)
        kl_loss = - 0.5 * K.sum(1 + z_log_var_encoded -
                                K.square(z_mean_encoded) -
                                K.exp(z_log_var_encoded), axis=-1)
        return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))

    def call(self, inputs):
        x = inputs[0]
        x_decoded = inputs[1]
        loss = self.vae_loss(x, x_decoded)
        self.add_loss(loss, inputs=inputs)
        # We won't actually use the output.
        return x


class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa

    # Behavior on each epoch
    def on_epoch_end(self, epoch, logs={}):
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)



# Process data

# Split 10% test set randomly
test_set_percent = 0.1
rnaseq_test_df = rnaseq_df.sample(frac=test_set_percent)
rnaseq_train_df = rnaseq_df.drop(rnaseq_test_df.index)

# Input place holder for RNAseq data with specific input size
rnaseq_input = Input(shape=(original_dim, ))

# ~~~~~~~~~~~~~~~~~~~~~~
# ENCODER
# ~~~~~~~~~~~~~~~~~~~~~~
# Depending on the depth of the model, the input is eventually compressed into
# a mean and log variance vector of prespecified size. Each layer is
# initialized with glorot uniform weights and each step (dense connections,
# batch norm,and relu activation) are funneled separately
#
# Each vector of length `latent_dim` are connected to the rnaseq input tensor
# In the case of a depth 2 architecture, input_dim -> latent_dim -> latent_dim2

if depth == 1:
    z_shape = latent_dim
    z_mean_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(rnaseq_input)
    z_log_var_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(rnaseq_input)

elif depth == 2:
    z_shape = latent_dim
    hidden_dense_linear = Dense(hidden_dim, kernel_initializer='glorot_uniform')(rnaseq_input)
    hidden_dense_batchnorm = BatchNormalization()(hidden_dense_linear)
    hidden_encoded = Activation('relu')(hidden_dense_batchnorm)

    z_mean_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(hidden_encoded)
    z_log_var_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(hidden_encoded)

elif depth == 3:
    z_shape = latent_dim
    hidden_dense_linear = Dense(hidden_dim, kernel_initializer='glorot_uniform')(rnaseq_input)
    hidden_dense_batchnorm = BatchNormalization()(hidden_dense_linear)
    hidden_encoded = Activation('relu')(hidden_dense_batchnorm)

    hidden2_dense_linear = Dense(hidden_dim_2, kernel_initializer='glorot_uniform')(hidden_encoded)
    hidden2_dense_batchnorm = BatchNormalization()(hidden2_dense_linear)
    hidden2_encoded = Activation('relu')(hidden2_dense_batchnorm)

    z_mean_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(hidden2_encoded)
    z_log_var_dense_linear = Dense(latent_dim, kernel_initializer='glorot_uniform')(hidden2_encoded)


z_mean_dense_batchnorm = BatchNormalization()(z_mean_dense_linear)
z_mean_encoded = Activation('relu')(z_mean_dense_batchnorm)

z_log_var_dense_batchnorm = BatchNormalization()(z_log_var_dense_linear)
z_log_var_encoded = Activation('relu')(z_log_var_dense_batchnorm)

# return the encoded and randomly sampled z vector
# Takes two keras layers as input to the custom sampling function layer with a
# latent_dim` output
z = Lambda(sampling,
           output_shape=(z_shape, ))([z_mean_encoded, z_log_var_encoded])

# ~~~~~~~~~~~~~~~~~~~~~~
# DECODER
# ~~~~~~~~~~~~~~~~~~~~~~
# The layers are different depending on the prespecified depth.
#
# Single layer: glorot uniform initialized and sigmoid activation.
# Double layer: relu activated hidden layer followed by sigmoid reconstruction
if depth == 1:
    decoder_to_reconstruct = Dense(original_dim,
                                   kernel_initializer='glorot_uniform',
                                   activation='sigmoid')
elif depth == 2:
    decoder_to_reconstruct = Sequential()
    decoder_to_reconstruct.add(Dense(hidden_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='relu',
                                     input_dim=latent_dim))
    decoder_to_reconstruct.add(Dense(original_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='sigmoid'))
elif depth == 3:
    decoder_to_reconstruct = Sequential()
    decoder_to_reconstruct.add(Dense(hidden_dim_2, 
                                     kernel_initializer='glorot_uniform',
                                     activation='relu',
                                     input_dim=latent_dim))
    decoder_to_reconstruct.add(Dense(hidden_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='relu',
                                     input_dim=hidden_dim_2))
    decoder_to_reconstruct.add(Dense(original_dim,
                                     kernel_initializer='glorot_uniform',
                                     activation='sigmoid'))

rnaseq_reconstruct = decoder_to_reconstruct(z)

# ~~~~~~~~~~~~~~~~~~~~~~
# CONNECTIONS
# ~~~~~~~~~~~~~~~~~~~~~~
adam = optimizers.Adam(lr=learning_rate)
vae_layer = CustomVariationalLayer()([rnaseq_input, rnaseq_reconstruct])
vae = Model(rnaseq_input, vae_layer)
vae.compile(optimizer=adam, loss=None, loss_weights=[beta])

# fit Model
hist = vae.fit(np.array(rnaseq_train_df),
               shuffle=True,
               epochs=epochs,
               batch_size=batch_size,
               validation_data=(np.array(rnaseq_test_df),
                                np.array(rnaseq_test_df)),
               callbacks=[WarmUpCallback(beta, kappa)])

# Save models
encoder = Model(rnaseq_input, z_mean_encoded)
encoded_rnaseq_df = encoder.predict_on_batch(rnaseq_df)
encoded_rnaseq_df = pd.DataFrame(encoded_rnaseq_df, index=rnaseq_df.index)

encoded_rnaseq_df.columns.name = 'sample_id'
encoded_rnaseq_df.columns = encoded_rnaseq_df.columns + 1
encoded_file = os.path.join('models', 
                            input_file + '_depth' + str(depth) + '_dm' + str(latent_dim) + '_onehidden_warmup_batchnorm.tsv')
encoded_rnaseq_df.to_csv(encoded_file, sep='\t')

# build a generator that can sample from the learned distribution
decoder_input = Input(shape=(latent_dim, ))  # can generate from any sampled z vector
_x_decoded_mean = decoder_to_reconstruct(decoder_input)
decoder = Model(decoder_input, _x_decoded_mean)

encoder_model_file = os.path.join('models', 
                                  input_file + '_depth' + str(depth) + '_dm' + str(latent_dim) + '_encoder_onehidden_vae.hdf5')
decoder_model_file = os.path.join('models', 
                                  input_file + '_depth' + str(depth) + '_dm' + str(latent_dim) + '_decoder_onehidden_vae.hdf5')

encoder.save(encoder_model_file)
decoder.save(decoder_model_file)

# Perform t-SNE on VAE encoded_features
tsne = manifold.TSNE(n_components=2, init='pca', random_state=0, perplexity=20,
                     learning_rate=300, n_iter=400)
tsne_out = tsne.fit_transform(encoded_rnaseq_df)
tsne_out = pd.DataFrame(tsne_out, columns=['1', '2'])
tsne_out.index = encoded_rnaseq_df.index
tsne_out.index.name = 'id'
tsne_out_file = os.path.join('features', 
                             input_file + '_depth' + str(depth) + '_dm' + str(latent_dim) + '_tybalt_tsne_features.tsv')
tsne_out.to_csv(tsne_out_file, sep='\t')

# Perform t-SNE on single cell gene expression features
if tsne_RNAseq == "True":
    tsne = manifold.TSNE(n_components=2, init='pca', random_state=0, perplexity=20,
                     learning_rate=300, n_iter=400)
    tsne_out = tsne.fit_transform(rnaseq_df)
    tsne_out = pd.DataFrame(tsne_out, columns=['1', '2'])
    tsne_out.index = rnaseq_df.index
    tsne_out.index.name = 'id'
    tsne_out_file = os.path.join('features', 
                             input_file + '_depth' + str(depth) + '_rnaseq_tsne_features.tsv')
    tsne_out.to_csv(tsne_out_file, sep='\t')


# Save training performance
output_filename = os.path.join('models', input_file + '_depth' + str(depth) + '_dm' + str(latent_dim) + '_training.perf.tsv')
history_df = pd.DataFrame(hist.history)
history_df = history_df.assign(learning_rate=learning_rate)
history_df = history_df.assign(batch_size=batch_size)
history_df = history_df.assign(epochs=epochs)
history_df = history_df.assign(kappa=kappa)
history_df = history_df.assign(seed=seed)
history_df = history_df.assign(depth=depth)
history_df = history_df.assign(first_layer=first_layer)
history_df = history_df.assign(latent_dim=latent_dim)
history_df.to_csv(output_filename, sep='\t')


