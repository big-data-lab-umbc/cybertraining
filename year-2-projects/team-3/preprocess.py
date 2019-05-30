import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from glob import glob
from os.path import join, expanduser
from sklearn.preprocessing import StandardScaler
from ipywidgets import interact
import ipywidgets as widgets

######## preprocess data

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
print('Shape of train_out',train_out.shape)
indices=np.where(train_out>out_threshold)[0]
train_out_pos=train_out[indices]
train_out_new=np.tile(train_out_pos,18)
train_out_new=np.concatenate((train_out,train_out_new),axis=0)
print(train_out_pos)
print(type(train_data))
train_data_pos=train_data.loc[indices,:]
#print(train_data_pos)
#print('Indices',indices)
test_data = csv_data.loc[valid_dates >= train_test_date, input_columns]
test_out = np.where(csv_data.loc[valid_dates >= train_test_date, output_column] > out_threshold, 1, 0)
print('Shape of test_out',test_out.shape)
print(f"Training strong rotation examples: {100 * train_out_new.sum() / train_out_new.size: 0.3f}%")
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
#    print(ds)
    # Stack the variables in the order listed within a given file
    in_data.append(np.stack([ds[v].values for v in in_vars], axis=-1))
    out_data.append(np.stack([ds[v].values for v in out_vars], axis=-1))
    
#    print(out_data)
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
#train_norm_2d_pos, scaling_values=normalize_multivariate_data(all_in_data[train_out>out_threshold])
train_norm_2d, scaling_values = normalize_multivariate_data(all_in_data[valid_dates < train_test_date])
#print(type(train_norm_2d))
#print('STOP HERE')
train_norm_2d_pos=train_norm_2d[indices,:,:,:]
train_norm_2d_new=np.tile(train_norm_2d_pos,(18,1,1,1))
print(train_norm_2d.shape)
print(train_norm_2d_new.shape)
print('STOP HERE')
train_norm_2d_new=np.concatenate((train_norm_2d,train_norm_2d_new),axis=0)

#print(train_norm_2d.shape)
#print(train_norm_2d_new.shape)
print('STOP HERE')
test_norm_2d, _ = normalize_multivariate_data(all_in_data[valid_dates >= train_test_date], scaling_values=scaling_values)
print(scaling_values)


np.save('data/train_out.npy', train_out)
np.save('data/test_out.npy', test_out)
np.save('data/train_norm_2d.npy', train_norm_2d)
np.save('data/test_norm_2d.npy', test_norm_2d)
np.save('data/train_out_new.npy', train_out_new)
np.save('data/train_norm_2d_new.npy', train_norm_2d_new)
