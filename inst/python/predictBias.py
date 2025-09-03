import os
import sys
import h5py
import pandas as pd
import numpy as np
import tqdm
import pickle
import multiprocessing as mp
import scipy.stats as ss
from datetime import datetime
from keras.models import load_model
from keras.models import Sequential
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
    
###############################################
# Use the model to predict Tn5 bias for regions #
###############################################
def main():
    args = pd.read_csv("args.txt", header = None)
    main_dir = args.values[0][0]
    data_dir = args.values[1][0]
    model_use = args.values[2][0]
    chunk_size = int(args.values[3][0])

    context_radius = 50

    # Load Tn5 bias model
    model = load_model(main_dir + model_use, compile=False)
    
    # Load sequences of regions
    print("Loading sequences of regions")
    region_seqs = pd.read_csv(data_dir + "regionSeqs.txt", header = None)
    region_seqs = np.transpose(region_seqs.values)[0]
    
    # Specify the radius of sequence context and the width of the region
    context_len = 2 * context_radius + 1
    region_width = len(region_seqs[0]) - 2 * context_radius
    
    # To reduce memory usage, we chunk the region list in to smaller chunks
    starts = np.arange(len(region_seqs), step = chunk_size)
    if len(starts) > 1:
        starts = starts[:(len(starts) - 1)]
        ends = starts + chunk_size
        ends[len(ends) - 1] = len(region_seqs) + 1
    else:
        starts = [0]
        ends = [len(region_seqs) + 1]
    
    # Create folder for storing intermediate results for each chunk
    os.system("mkdir " + data_dir + "chunked_bias_pred")
    
    # Go through all chunks and predict Tn5 bias
    print("Predicting Tn5 bias for regions")
    for i in tqdm.tqdm(range(len(starts))):
    
        if os.path.exists(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl"):
            continue
    
        print("Processing chunk No." + str(i) + " " + datetime.now().strftime("%H:%M:%S"))
        chunk_seqs = region_seqs[starts[i]:ends[i]]
        
        # Encode sequences in the current chunk into matrices using one-hot encoding
        print("Encoding sequence contexts")
        with mp.Pool(2) as pool:
           chunk_onehot = list(tqdm.tqdm(pool.imap(region_onehot_encode, chunk_seqs), total = len(chunk_seqs)))
        
        # Use neural network to predict bias
        # For sequences containing "N", the encoded matrix will be a empty zero matrix,
        # in such cases we assign 1s as bias values
        print("Predicting bias")
        pred_bias = np.array([np.transpose(model.predict(chunk_onehot[j]))[0] \
                              if (np.sum(chunk_onehot[j]) > 1) else np.ones(region_width) \
                              for j in tqdm.tqdm(range(len(chunk_onehot)), total = chunk_size)], 
                            dtype=object)
        
        # Reverse transform the predicted values to the original scale
        pred_bias = np.power(10, (pred_bias - 0.5) * 2) - 0.01
    
        # Save intermediate result
        with open(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl", "wb") as f:
            pickle.dump(pred_bias, f)
    
    # Integrate results from all chunks
    pred_bias_all = None
    for i in tqdm.tqdm(range(len(starts))):
    
        # Load intermediate result
        with open(data_dir + "chunked_bias_pred/chunk_" + str(i) + ".pkl", "rb") as f:
            pred_bias = pickle.load(f)
    
        # Reverse transform the predicted values to the original scale
        if pred_bias_all is None:
            pred_bias_all = pred_bias
        else:
            pred_bias_all = np.concatenate([pred_bias_all, pred_bias])
    
    # Write results of the current batch to file
    # output_path = data_dir + "predBias.h5"
    # if os.path.isfile(output_path):
    #     os.system("rm " + output_path)
    # hf = h5py.File(output_path, 'w')
    # hf.create_dataset("predBias", data = pred_bias_all)
    # hf.close()

    max_len = max(len(row) for row in pred_bias_all)
    pred_bias_all = [np.pad(row, (0, max_len - len(row)), constant_values=np.nan) for row in pred_bias_all]
    pred_bias_all = np.array(pred_bias_all, dtype=object)

    np.savetxt(data_dir + "pred_bias.txt", pred_bias_all)

    # Remove intermediate files
    os.system("rm -r " + data_dir + "chunked_bias_pred")


if __name__ == "__main__":
    import multiprocessing as _mp
    try:
        _mp.set_start_method("spawn")
    except RuntimeError:
        pass
    main()
    
