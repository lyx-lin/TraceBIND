import os
import sys
import h5py
import pandas as pd
import numpy as np
import tqdm
import random
import pickle
import multiprocessing as mp
import scipy.stats as ss
import tensorflow as tf
from tensorflow.keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from datetime import datetime
from keras.models import load_model
from keras.models import Sequential
from keras.optimizers import Adam
from keras.layers import *
from keras.models import Model
from keras import backend as K

print("Python CWD:", os.getcwd())

def onehot_encode(seq):
    mapping = pd.Series(index = ["A", "C", "G", "T"], data = [0, 1, 2, 3])
    bases = [base for base in seq]
    base_inds = mapping[bases]
    onehot = np.zeros((len(bases), 4))
    onehot[np.arange(len(bases)), base_inds] = 1
    return onehot

def region_onehot_encode(region_seq, context_radius = 50):
    
    # Calculate length of every local sub-sequence
    context_len = 2 * context_radius + 1
    
    # If region width is L, then region_seq should be a string of length L + 2 * context_len
    region_width = len(region_seq) - 2 * context_radius
    
    if "N" in region_seq:
        return np.zeros((region_width, 4))
    else:
        # First encode the whole region sequence 
        # This prevents repetitive computing for overlapping sub-sequences
        region_onehot = np.array(onehot_encode(region_seq))
        
        # Retrieve encoded sub-sequences by subsetting the larger encoded matrix
        region_onehot = np.array([region_onehot[i : (i + context_len), :] for i in range(region_width)])
        
        return region_onehot



def main():
    args = pd.read_csv("args.txt", header = None)
    obsbias_path = args.values[0][0]
    PRINT_model_path = args.values[1][0]
    finetuned_model_save_path = args.values[2][0]
    finetuned_model_name = args.values[3][0]
    n_jobs = int(args.values[4][0])

    context_radius = 50
    chunk_size = 2000

    bias_data = pd.read_csv(obsbias_path, sep = "\t")
    print("One-hot encoding of sequence contexts")
    seqs = bias_data.loc[:, "context"]
    with mp.get_context("spawn").Pool(n_jobs) as pool:
        onehot_seqs = list(tqdm.tqdm(pool.imap(onehot_encode, seqs), total=len(seqs)))
    onehot_seqs = np.array(onehot_seqs)

    target = bias_data.loc[:, "obs_bias"].values
    target = np.log10(target + 0.01)
    target = (target / 2) + 0.5

    mapped_regions = bias_data.loc[:, "BACInd"].values
    n_regions = len(np.unique(mapped_regions))

    np.random.seed(1307)
    region_inds = np.unique(mapped_regions)
    np.random.shuffle(region_inds)
    training_region_inds = region_inds[:int(n_regions * 0.95)]
    test_region_inds = region_inds[int(n_regions * 0.95):]
    print(test_region_inds)

    training_inds = np.array([i for i in range(len(mapped_regions)) if \
                    mapped_regions[i] in training_region_inds])
    test_inds = np.array([i for i in range(len(mapped_regions)) if \
                    mapped_regions[i] in test_region_inds])

    training_data = onehot_seqs[training_inds]
    training_target = target[training_inds]
    test_data = onehot_seqs[test_inds]
    test_target = target[test_inds]

    np.random.shuffle(training_inds)
    training_data = onehot_seqs[training_inds]
    training_target = target[training_inds]

    training_data = training_data.reshape((-1, 101, 4))
    training_target = training_target.reshape((-1, 1))

    print("Training Tn5 bias model")
    inputs = Input(shape = (np.shape(test_data[0])), name='input_1')
    conv_1 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1, name='conv1d')(inputs)
    maxpool_1 = MaxPooling1D(name='max_pooling1d')(conv_1)
    conv_2 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1, name='conv1d_1')(maxpool_1)
    maxpool_2 = MaxPooling1D(name='max_pooling1d_1')(conv_2)
    conv_3 = Conv1D(32, 5, padding = 'same', activation = 'relu', strides = 1, name='conv1d_2')(maxpool_2)
    maxpool_3 = MaxPooling1D(name='max_pooling1d_2')(conv_3)
    flat = Flatten(name='flatten')(maxpool_3)
    fc = Dense(32, activation = "relu", name='dense')(flat)
    # fc = BatchNormalization()(fc)
    # fc = Dropout(0.1)(fc) 
    # fc = Dense(16, activation = "relu", name='dense_2')(flat)
    # fc = BatchNormalization()(fc)
    # fc = Dropout(0.2)(fc) 
    out = Dense(1, activation = "linear", name='dense_1')(fc)
    model = Model(inputs=inputs, outputs=out) 
    # model.layers[1].trainable = False
    model.layers[3].trainable = False
    # model.layers[5].trainable = False
    # model.layers[8].trainable = False

    model.summary()
    optimizer = Adam(learning_rate=0.0001)

    model.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['mse'])
    model.load_weights(PRINT_model_path, by_name=True, skip_mismatch=True
                      )

    test_pred = np.transpose(model.predict(test_data))[0]
    test_pred_rev = np.power(10.0, (test_pred - 0.5) * 2) - 0.01
    test_target_rev = np.power(10.0, (test_target - 0.5) * 2) - 0.01
    # print("Pearson correlation = " + str(ss.pearsonr(test_target, test_pred)[0]))
    # print("Pearson correlation = " + str(ss.pearsonr(test_target_rev, test_pred_rev)[0]))
    PRINT_cor = str(ss.pearsonr(test_target, test_pred)[0])
    callbacks = []
    lr_cb = ReduceLROnPlateau(monitor='val_loss', patience=4, factor=0.1, verbose=True)
    callbacks.append(lr_cb)
    es_cb = EarlyStopping(monitor='val_loss', patience=6, verbose=True)
    callbacks.append(es_cb)
    model.fit(training_data, 
              training_target, 
              batch_size=64,
              epochs = 300, 
              shuffle=True,
              callbacks=callbacks,
              validation_split=0.1,
              verbose=0)   

    # Model evaluation on the test set
    print("Evaluating performance on the test set")
    # plt_ind = np.random.choice(np.arange(len(onehot_seqs)), 10000)
    print(PRINT_cor)
    test_pred = np.transpose(model.predict(test_data))[0]
    test_target_rev = np.power(10, (test_target - 0.5) * 2) - 0.01
    test_pred_rev = np.power(10, (test_pred - 0.5) * 2) - 0.01

    print("Pearson correlation before finetuning = " + PRINT_cor, flush=True)

    print("Pearson correlation after finetuning = " + str(ss.pearsonr(test_target, test_pred)[0]), flush=True)

    model.save(finetuned_model_save_path + finetuned_model_name) 



if __name__ == "__main__":
    import multiprocessing as _mp
    try:
        _mp.set_start_method("spawn")
    except RuntimeError:
        pass
    main()
