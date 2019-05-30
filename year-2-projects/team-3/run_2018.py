
# AUTHORS: CHARLIE BECKER, BIN WANG, WILL MAYFIELD, SARAH MURPHY
# DATE: 05/22/2019

# THIS SCRIPT DEMONSTRATES A METHOD OF TUNING HYPERPARAMETERS OF A DEEP 
# NEURAL NETWORK IN PARALLEL IN AN HPC ENVIRONMENT USING A COMBINATION 
# OF POPULAR PYTHON MODULES - DASK, SCIKIT-LEARN AND KERAS. DATA AND BASE
# CONVOLUTIONAL MODEL STRUCTURE IS BORROWED FROM THE 2019 AMS SHORT COURSE
# ON MAHCINE LEARNING IN PYTHON TAUGHT BY JOHN GAGNE FROM NCAR. THE SHORT 
# COURSE GITHUB CAN BE FOUND AT https://github.com/djgagne/ams-ml-python-course

from sklearn.model_selection import GridSearchCV
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
from keras.preprocessing.image import ImageDataGenerator
from dask_jobqueue import SLURMCluster
from dask.distributed import Client
from dask.distributed import progress
from keras.wrappers.scikit_learn import KerasClassifier, KerasRegressor
#from imblearn.over_sampling import RandomOverSampler

####################### INITIATE CLUSTER ##################################

cluster = SLURMCluster(cores=1, memory='48 GB', job_extra=['--account=pi_gobbert','--partition=high_mem','--qos=normal+','--time=03:00:00','--exclusive'])
cluster.scale(5)

client = Client(cluster)
print(cluster.job_script())
####################### OUTER HYPERPARAMTER ###############################

bs = [256,512,1024,2048,4096]
start = time.time()

###################### RETRIEVE MODEL METRICS #############################

def calc_verification_scores(test_labels,predictions):
    
    model_auc = roc_auc_score(test_labels, predictions)
    model_brier_score = mean_squared_error(test_labels, predictions)
    climo_brier_score = mean_squared_error(test_labels, np.ones(test_labels.size) * test_labels.sum() / test_labels.size)
    model_brier_skill_score = 1 - model_brier_score / climo_brier_score
    print(f"AUC: {model_auc:0.3f}")
    print(f"Brier Score: {model_brier_score:0.3f}")
    print(f"Brier Score (Climatology): {climo_brier_score:0.3f}")
    print(f"Brier Skill Score: {model_brier_skill_score:0.3f}")
    return model_auc

###################### BUILD MODEL STRUCTURE ##############################

def create_model(learning_rate='0.001'):
    
    # Deep convolutional neural network
    num_conv_filters = 8
    filter_width = 5
    conv_activation = "relu"
    
    # Input data in shape (instance, y, x, variable)
    conv_net_in = Input(shape=(32, 32, 3))  # train_norm_2d.shape[1:]
    
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
    
    # Compile
    conv_model.compile(opt, "binary_crossentropy", metrics=['accuracy'])

    print(conv_model.summary())
    return conv_model

#################### LOAD DATA AND DISTRIBUTE MODEL ########################

def train_model(param):
    
    train_out = np.load('data/train_out.npy')
    test_out = np.load('data/test_out.npy')
    train_norm_2d = np.load('data/train_norm_2d.npy')
    test_norm_2d = np.load('data/test_norm_2d.npy')
    train_norm_2d_new = np.load('data/train_norm_2d_new.npy')
    train_out_new = np.load('data/train_out_new.npy')
    out_threshold = 0.005
    
    ####################### GRID HYPERPARAMETERS ###########################
    
    KFolds = 3
    e = [5,10,15]
    learning_rate = [0.01,0.005,0.001]
    param_grid = dict(epochs=e,learning_rate=learning_rate)

    ####################### CV GRID SEARCH #################################
    augmentation=False

    t0 = time.time()
    l = []
    
    model = KerasClassifier(build_fn=create_model, batch_size = param)
    
    if augmentation==True:
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
       print("Running augmented training now, with augmentation")
       history=model.fit_generator(datagen.flow(train_norm_2d_new, train_out_new, batch_size=param,shuffle=True),steps_per_epoch=len(train_norm_2d_new)/param,epochs=10)
       indices_test=np.where(test_out>out_threshold)[0]
       test_out_pos=test_out[indices_test]
       test_out_new=np.tile(test_out_pos,18)
       test_out_new=np.concatenate((test_out,test_out_new),axis=0)
       test_norm_2d_pos=test_norm_2d[indices_test,:,:,:]
       test_norm_2d_new=np.tile(test_norm_2d_pos,(18,1,1,1))
       test_norm_2d_new=np.concatenate((test_norm_2d,test_norm_2d_new),axis=0)
       t1 = time.time()
       return calc_verification_scores(model, test_norm_2d_new, test_out_new), t1-t0


    else:
      print("Running regular training, no augmentation")    
      #ros=RandomOverSampler(random_state=0)
      #train_norm_2d_resampled, train_out_resampled=ros.fit_resample(train_norm_2d,train_out)
      
      start_time = time.time()
      grid = GridSearchCV(estimator=model, param_grid = param_grid, cv = KFolds,  n_jobs = -1)
      grid_result = grid.fit(train_norm_2d_new, train_out_new)
      
      means = grid_result.cv_results_['mean_test_score']
      stds = grid_result.cv_results_['std_test_score']
      params = grid_result.cv_results_['params']  
      
      end_time = time.time()
      run_time = end_time - start_time
      #preds = model.predict(test_norm_2d)
      
      for mean, stdev, p in zip(means, stds, params):
        x = [("%f (%f) with: %r, batch: %r ,time: %r" % (mean, stdev, p, param, run_time))]
        l.append(x)
      
      for item in l:
        print(item)
      #calc_verification_scores(test_out,preds)
      return l

################### DISTRIBUTE MODELS AND GATHER RESULTS ##################

models = client.map(train_model, bs)       
results = client.gather(models)
for item in results:
    print(*item, sep='\n')        
end=time.time()
print('Total time used is ', np.round(end-start,3))
client.close()
