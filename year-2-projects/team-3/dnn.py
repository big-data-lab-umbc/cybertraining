import time
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from glob import glob
from os.path import join, expanduser
from sklearn.preprocessing import StandardScaler
from ipywidgets import interact
import ipywidgets as widgets
from keras.models import Model, save_model, load_model
from keras.layers import Dense, Activation, Conv2D, Input, AveragePooling2D, MaxPooling2D, Flatten, LeakyReLU, Dropout
from keras.layers import SpatialDropout2D
from keras.optimizers import SGD, Adam
from keras.regularizers import l2
import keras.backend as K
from scipy.ndimage import gaussian_filter
from sklearn.metrics import mean_squared_error, roc_auc_score

csv_path = "data/track_data_ncar_ams_3km_csv_small/"
csv_files = sorted(glob(join(csv_path, "*.csv")))
csv_data_list = []
for csv_file in csv_files:
    print(csv_file)
    csv_data_list.append(pd.read_csv(csv_file))
csv_data = pd.concat(csv_data_list, ignore_index=True)
del csv_data_list

# input variables
input_columns = ["REFL_COM_mean", "U10_mean", "V10_mean", "T2_mean"]
output_column = "RVORT1_MAX-future_max"
# vorticity threshold in s-1
out_threshold = 0.005
train_test_date = pd.Timestamp("2015-01-01")
valid_dates = pd.DatetimeIndex(csv_data["Valid_Date"])
# Extract the input training data to the neural network
train_data = csv_data.loc[valid_dates < train_test_date, input_columns]
# Label strong rotation as 1 and weak or no rotation as 0
train_out = np.where(csv_data.loc[valid_dates < train_test_date, output_column] > out_threshold, 1, 0)
# 
test_data = csv_data.loc[valid_dates >= train_test_date, input_columns]
test_out = np.where(csv_data.loc[valid_dates >= train_test_date, output_column] > out_threshold, 1, 0)
print(f"Training strong rotation examples: {100 * train_out.sum() / train_out.size: 0.3f}%")
print(f"Training strong rotation example: {100 * test_out.sum() / test_out.size: 0.3f}%")

scaler = StandardScaler()
train_norm = scaler.fit_transform(train_data)
test_norm = scaler.transform(test_data)
for i, input_col in enumerate(input_columns):
    print(f"{input_columns[i]:13s} Mean: {scaler.mean_[i]:0.3f} SD: {scaler.scale_[i]:0.3f}")

# Create a list of storm netCDF4 files. If the files are not in your home directory,
# Change the start_path variable to the appropriate location
#start_path = join(expanduser("~"), "ams-ml-python-course")
storm_files = sorted(glob("data/track_data_ncar_ams_3km_nc_small/*.nc"))
print(storm_files[0])
run_times = []
valid_times = []
# List of input variables
in_vars = ["REFL_COM_curr",
           "U10_curr", "V10_curr"]
# List of output variables
out_vars = ["RVORT1_MAX_future"]
in_data = []
out_data = []
# Loop through each storm file and extract the relevant variables
for storm_file in storm_files:
    # Extract run time from the filename
    run_time = pd.Timestamp(storm_file.split("/")[-1].split("_")[1])
    # If you want to ignore certain run_dates, encapsulate the remaining lines in an if statement
    ds = xr.open_dataset(storm_file)
    # Stack the variables in the order listed within a given file
    in_data.append(np.stack([ds[v].values for v in in_vars], axis=-1))
    out_data.append(np.stack([ds[v].values for v in out_vars], axis=-1))
    # Extract the valid times
    valid_times.append(ds["time"].values)
    # Extract the run times and match each run time with each patch
    run_times.append([run_time] * in_data[-1].shape[0])
    print(run_time)
    ds.close()
# Stack the  data into single arrays instead of lists of arrays
all_in_data = np.vstack(in_data)
all_out_data = np.vstack(out_data)
all_run_times = np.concatenate(run_times)
all_valid_times = np.concatenate(valid_times)
# Deallocate the lists of arrays to save memory
del in_data[:], out_data[:], run_times[:], valid_times[:]
del in_data, out_data, run_times, valid_times


def normalize_multivariate_data(data, scaling_values=None):
    """
    Normalize each channel in the 4 dimensional data matrix independently.

    Args:
        data: 4-dimensional array with dimensions (example, y, x, channel/variable)
        scaling_values: pandas dataframe containing mean and std columns

    Returns:
        normalized data array, scaling_values
    """
    normed_data = np.zeros(data.shape, dtype=data.dtype)
    scale_cols = ["mean", "std"]
    if scaling_values is None:
        scaling_values = pd.DataFrame(np.zeros((data.shape[-1], len(scale_cols)), dtype=np.float32),
                                      columns=scale_cols)
    for i in range(data.shape[-1]):
        scaling_values.loc[i, ["mean", "std"]] = [data[:, :, :, i].mean(), data[:, :, :, i].std()]
        normed_data[:, :, :, i] = (data[:, :, :, i] - scaling_values.loc[i, "mean"]) / scaling_values.loc[i, "std"]
    return normed_data, scaling_values

train_norm_2d, scaling_values = normalize_multivariate_data(all_in_data[valid_dates < train_test_date])
test_norm_2d, _ = normalize_multivariate_data(all_in_data[valid_dates >= train_test_date], scaling_values=scaling_values)
print(scaling_values)

start=time.time()
# Train a simple 1-layer dense  neural network to predict the probability of strong rotation
# set the batch size 
batch_size = 128
# an epoch is one pass through the training data
num_epochs = 10
# number of neurons in the hidden layer
hidden_neurons = 10
# learning rate
learning_rate = 0.01
# Set the optimizer to SGD or change it to an Adam object
optimizer = SGD(lr=learning_rate, momentum=0.9)
# loss can be binary_crossentropy or mse
loss = "binary_crossentropy"
# Create an input layer
small_net_input = Input(shape=(len(input_columns),))
# Hidden layer with a hyperbolic tangent (tanh) activation function
small_net_hidden = Dense(hidden_neurons, activation="tanh")(small_net_input)
# Output layers with a sigmoid activation function
small_net_out = Dense(1, activation="sigmoid")(small_net_hidden)
# combine the layers into a model object
small_model = Model(small_net_input, small_net_out)
# compile the model object to instantiate all the tensorflow connections
small_model.compile(optimizer, loss=loss)
# The batch size should be divisible by the number of training examples. 
# We will subject the remainder to ensure this
batch_diff = train_norm.shape[0] % batch_size
# Create training indices for random batch sampling later
batch_indices = np.arange(train_norm.shape[0] - batch_diff)
num_batches = int(batch_indices.size / batch_size)
# Extract the weights from the network
# hw = hidden weights, hb = hidden bias, ow = output weights ,ob = output bias
hw, hb, ow, ob = small_model.get_weights()
# Store the values of the hidden and output weights after each batch
hw_series = np.zeros([num_batches * num_epochs] + list(hw.shape))
ow_series = np.zeros([num_batches * num_epochs] + list(ow.shape))
i = 0
# Each epoch is one pass through the training data
for e in range(num_epochs):
    print(e)
    # At the beginning of each epoch, shuffle the order of the training examples
    batch_indices = np.random.permutation(batch_indices)
    for b in range(num_batches):
        hw_series[i], hb, ow_series[i], ob = small_model.get_weights()
        small_model.train_on_batch(train_norm[batch_indices[b * batch_size: (b + 1) * batch_size]],
                                   train_out[batch_indices[b * batch_size: (b + 1) * batch_size]])
        
        i += 1

small_model.fit(train_norm, train_out, epochs=10, batch_size=128)

small_model.summary()

#calc_verification_scores(small_model, test_norm_2d, test_out)


small_model_preds = small_model.predict(test_norm).ravel()
small_model_auc = roc_auc_score(test_out, small_model_preds)
small_model_brier_score = mean_squared_error(test_out, small_model_preds)
climo_brier_score = mean_squared_error(test_out, np.ones(test_out.size) * test_out.sum() / test_out.size)
small_model_brier_skill_score = 1 - small_model_brier_score / climo_brier_score
print(f"AUC: {small_model_auc:0.3f}")
print(f"Brier Score: {small_model_brier_score:0.3f}")
print(f"Brier Score (Climatology): {climo_brier_score:0.3f}")
print(f"Brier Skill Score: {small_model_brier_skill_score:0.3f}")
end=time.time()
print('Time used is', end-start)




