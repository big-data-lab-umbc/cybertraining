import random
from random import shuffle
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
#from imblearn.over_sampling import RandomOverSampler
from keras.preprocessing.image import ImageDataGenerator
import skimage
from skimage import transform

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
indices_new_shuffled=np.arange(train_out_new.shape[0])
print('Print indices_new',indices_new_shuffled)
random.shuffle(indices_new_shuffled)
print('Print indices_new_shuffled',indices_new_shuffled)
train_out_new_shuffled=train_out_new[indices_new_shuffled]
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
train_norm_2d_new_shuffled=train_norm_2d_new[indices_new_shuffled,:,:,:]
train_norm_2d_new_shuffled_transformed=train_norm_2d_new_shuffled
for ii in np.arange(train_norm_2d_new.shape[0]):
    train_norm_2d_new_shuffled_transformed[ii,:,:,:]=skimage.transform.rotate(train_norm_2d_new_shuffled[ii,:,:,:],2,resize=False,center=None,order=1,mode='constant',preserve_range=True)
    print(train_norm_2d_new_shuffled_transformed[ii,:,:,:]) 

#print(train_norm_2d.shape)
#print(train_norm_2d_new.shape)
print('STOP HERE')
test_norm_2d, _ = normalize_multivariate_data(all_in_data[valid_dates >= train_test_date], scaling_values=scaling_values)
print(scaling_values)

indices_test=np.where(test_out>out_threshold)[0]
test_out_pos=test_out[indices_test]
test_out_new=np.tile(test_out_pos,18)
test_out_new=np.concatenate((test_out,test_out_new),axis=0)
test_norm_2d_pos=test_norm_2d[indices_test,:,:,:]
test_norm_2d_new=np.tile(test_norm_2d_pos,(18,1,1,1))
test_norm_2d_new=np.concatenate((test_norm_2d,test_norm_2d_new),axis=0)
perm=np.arange(test_out_new.shape[0])
print('STOP STOP STOP',perm.shape)
np.random.shuffle(perm)
test_norm_2d_new_shuffled=test_norm_2d_new[perm,:,:,:]
test_out_new_shuffled=test_out_new[perm]

def calc_verification_scores(model, test_data, test_labels):
    preds = model.predict(test_data, batch_size=1024)
    model_auc = roc_auc_score(test_labels, preds)
    model_brier_score = mean_squared_error(test_labels, preds)
    climo_brier_score = mean_squared_error(test_labels, np.ones(test_labels.size) * test_labels.sum() / test_labels.size)
    model_brier_skill_score = 1 - model_brier_score / climo_brier_score
    print(f"AUC: {model_auc:0.3f}")
    print(f"Brier Score: {model_brier_score:0.3f}")
    print(f"Brier Score (Climatology): {climo_brier_score:0.3f}")
    print(f"Brier Skill Score: {model_brier_skill_score:0.3f}")


start=time.time()
# Deep convolutional neural network
num_conv_filters = 8
filter_width = 5
conv_activation = "relu"
learning_rate = 0.001
# Input data in shape (instance, y, x, variable)
conv_net_in = Input(shape=train_norm_2d.shape[1:])
# First 2D convolution Layer
conv_net = Conv2D(num_conv_filters, (filter_width, filter_width), padding="same")(conv_net_in)
conv_net = Activation(conv_activation)(conv_net)
# Average pooling takes the mean in a 2x2 neighborhood to reduce the image size
conv_net = AveragePooling2D()(conv_net)
# Second set of convolution and pooling layers
conv_net = Conv2D(num_conv_filters * 2, (filter_width, filter_width), padding="same")(conv_net)
conv_net = Activation(conv_activation)(conv_net)
conv_net = AveragePooling2D()(conv_net)
# Third set of convolution and pooling layers
conv_net = Conv2D(num_conv_filters * 4, (filter_width, filter_width), padding="same")(conv_net)
conv_net = Activation(conv_activation)(conv_net)
conv_net = AveragePooling2D()(conv_net)
# Flatten the last convolutional layer into a long feature vector
conv_net = Flatten()(conv_net)
# Dense output layer, equivalent to a logistic regression on the last layer
conv_net = Dense(1)(conv_net)
conv_net = Activation("sigmoid")(conv_net)
conv_model = Model(conv_net_in, conv_net)
# Use the Adam optimizer with default parameters
opt = Adam(lr=learning_rate)
conv_model.compile(opt, "binary_crossentropy")

conv_model.summary()

augmentation=True
on_the_fly=False
if augmentation==True:
   if on_the_fly==True:
      datagen = ImageDataGenerator(
             rotation_range=5,
             width_shift_range=0,
             height_shift_range=0,
             shear_range=0,
             zoom_range=0,
             horizontal_flip=True,
             fill_mode='nearest')

      #datagen.fit(train_norm_2d_new)
      print(train_norm_2d_new.shape)
      print(train_out_new.shape)
      print("Running augmented training now, with augmentation on the fly")
      history=conv_model.fit_generator(datagen.flow(train_norm_2d_new_shuffled, train_out_new_shuffled, batch_size=512,shuffle=True),steps_per_epoch=len(train_norm_2d_new)/512,epochs=10)
   else:
      print("Running training with pre-augmented data")
      conv_model.fit(train_norm_2d_new_shuffled_transformed, train_out_new_shuffled,batch_size=512,epochs=10)

   calc_verification_scores(conv_model, test_norm_2d_new_shuffled, test_out_new_shuffled)
else:
   print("Running regular training, no augmentation")    
   #ros=RandomOverSampler(random_state=0)
   #train_norm_2d_resampled, train_out_resampled=ros.fit_resample(train_norm_2d,train_out)
   conv_model.fit(train_norm_2d, train_out, batch_size=512, epochs=10)
   calc_verification_scores(conv_model, test_norm_2d, test_out)

end=time.time()
print('Time used is ', end-start)



