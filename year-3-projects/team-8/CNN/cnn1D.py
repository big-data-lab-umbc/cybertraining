import numpy as np
from readRUCsounding import readRUCsounding
from convert_variables import convert_variables
from interpolateRUCsounding import interpolate_height
from normalize_variables import normalize_data
from evaluate_model import evaluate_model
from imbalanced_data import *
from sklearn.model_selection import train_test_split
from tensorflow.keras.utils import to_categorical

## path to data
path='/umbc/xfs1/cybertrn/cybertraining2020/team8/research/RUCsoundings/files/'

## read data
RUCdata, RUCsounding = readRUCsounding(path)

## interpolate sounding data to fixed height grids
RUCsounding = interpolate_height(RUCsounding)

## convert tc,tdc to theta,q
#RUCsounding = convert_variables(RUCsounding)

## remove prs, z: after interpolation, these are implicit in other variables
RUCsounding = RUCsounding[:,:,[1, 2, 3, 4, 5]] 

## create label array
labels = np.array(RUCdata['eventtype']).reshape(len(RUCdata['eventtype']),1)
labels = to_categorical(labels)

## split soundings into training and testing data
training_data, testing_data, training_labels, testing_labels = \
    train_test_split(RUCsounding, labels,
                     test_size=0.2, random_state=0)

## normalize sounding data
training_data, testing_data = normalize_data(training_data, testing_data)

## refit training data to deal with imbalanced classes
training_data, training_labels = oversample_old(training_data, training_labels)

n_levels, n_variables, n_outputs = training_data.shape[1], \
                                   training_data.shape[2], \
                                   training_labels.shape[1]

## evaluate CNN for ten iterations
scores = list()
for r in range(10):
        _, score = evaluate_model(training_data, training_labels,
                                  testing_data,  testing_labels)
        score = score * 100.0
        print('>#%d: %.3f' % (r+1, score))
        scores.append(score)

m, s = np.mean(scores), np.std(scores)
print('Accuracy: %.3f%% (+/-%.3f)' % (m, s))
