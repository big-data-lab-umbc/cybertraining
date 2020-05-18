import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tensorflow import keras
from tensorflow.keras import models
from tensorflow.keras.callbacks import CSVLogger
from tensorflow.keras.callbacks import ModelCheckpoint
from create_model import create_model_v2 as create_model
from sklearn.metrics import confusion_matrix

def evaluate_model(training_data, training_labels, testing_data, testing_labels):

    ## evaluate model with: two 1D convolution layers
    #                       one Dropout layer,
    #                       one 1D max pooling layer, and
    #                       two Dense layers
   
    n_levels, n_variables, n_outputs = training_data.shape[1], \
                                       training_data.shape[2], \
                                       training_labels.shape[1]

    ## create CNN
    model = create_model(n_levels, n_variables, n_outputs)

    ## fit model
    filepath="weights-improvement-{epoch:02d}-{val_loss:.2f}.hdf5"
    checkpoint = ModelCheckpoint(filepath,monitor='val_loss', period=8,verbose=1)
    callbacks_list = [checkpoint]

    verbose, epochs, batch_size = 0, 50, 10
    history = model.fit(training_data, training_labels, validation_split=0.2, 
                        epochs=epochs, batch_size=batch_size, verbose=verbose,
                        shuffle=True, callbacks=callbacks_list)

    acc = history.history['accuracy']
    val_acc = history.history['val_accuracy']
    loss = history.history['loss']
    val_loss = history.history['val_loss']
    epochs = range(1, len(acc) + 1)
    plt.subplot(122)
    plt.plot(epochs, loss, 'b:', label='Training loss')
    plt.plot(epochs, val_loss, 'b', label='Validation loss')
    plt.plot(epochs, acc, 'r:', label='Training accuracy')
    plt.plot(epochs, val_acc, 'r', label='Validation accuracy')
    plt.ylim([0,2.0])
    plt.title('Training and validation accuracy/loss')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy/Loss')
    plt.legend()
    plt.savefig('image.png', bbox_inches='tight',dpi=300)
    plt.close()

    ## evaluate model
    (loss, accuracy) = model.evaluate(testing_data, testing_labels,
                                      batch_size=batch_size, verbose=verbose)

    preds = np.round(model.predict(testing_data),0)

    categorical_test_labels = pd.DataFrame(testing_labels).idxmax(axis=1)
    categorical_preds = pd.DataFrame(preds).idxmax(axis=1)

    CM = confusion_matrix(categorical_test_labels,
                          categorical_preds)
    print('Confusion Matrix: \n', CM)

    return loss, accuracy
