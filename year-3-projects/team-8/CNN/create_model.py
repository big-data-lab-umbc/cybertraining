from tensorflow import keras
from tensorflow.keras import models
from tensorflow.keras.layers import (Conv1D, Dropout, MaxPool1D, 
                                     Flatten, Dense, BatchNormalization, 
                                     LeakyReLU)
from tensorflow.keras.metrics import FalseNegatives

def create_model_v2(n_levels, n_variables, n_outputs):

    ## build model with: three 1D convolution layers,
    #                    three BN layers,
    #                    three 1D max pooling layers, and
    #                    two Dense layers

    model=models.Sequential()
    model.add(Conv1D(filters=32, kernel_size=5, 
                     activation=LeakyReLU(alpha=0.1), 
                     input_shape=(n_levels,
                                  n_variables)))
    model.add(BatchNormalization())
    model.add(MaxPool1D(pool_size=2))
    model.add(Conv1D(filters=64, kernel_size=3, 
                     activation=LeakyReLU(alpha=0.1)))
    model.add(BatchNormalization())
    model.add(MaxPool1D(pool_size=2))
    model.add(Conv1D(filters=128, kernel_size=3, 
                     activation=LeakyReLU(alpha=0.1)))
    model.add(BatchNormalization())
    model.add(MaxPool1D(pool_size=2))
#   model.add(Conv1D(filters=256, kernel_size=3, 
#                    activation=LeakyReLU(alpha=0.1)))
#   model.add(BatchNormalization())
#   model.add(MaxPool1D(pool_size=2))
    model.add(Flatten())
#   model.add(Dense(units=100, activation=LeakyReLU(alpha=0.1)))
    model.add(Dense(units=n_outputs, activation='softmax'))

    ## compile model
    model.compile(optimizer='sgd',
                    loss='mse',
#                   metrics=[FalseNegatives()])
                    metrics=['accuracy'])

    return model

def create_model_v1(n_levels, n_variables, n_outputs):

    ## build model with: two 1D convolution layers
    #                    one Dropout layer,
    #                    one 1D max pooling layer, and
    #                    two Dense layers

    model=models.Sequential()
    model.add(Conv1D(filters=64, kernel_size=3,
                     activation='relu',
                     input_shape=(n_levels,
                                  n_variables)))
    model.add(Conv1D(filters=64, kernel_size=3,
                     activation='relu'))
#   model.add(Dropout(0.3))
    model.add(MaxPool1D(pool_size=2))
    model.add(Flatten())
    model.add(Dense(units=100, activation='relu'))
    model.add(Dense(units=n_outputs, activation='softmax'))

    ## compile model
    model.compile(optimizer='adam',
                    loss='categorical_crossentropy',
                    metrics=['accuracy'])

#   print(model.summary())
    return model
