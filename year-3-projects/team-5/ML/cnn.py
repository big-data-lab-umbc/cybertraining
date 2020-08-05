import tensorflow.keras
from tensorflow.keras.datasets import fashion_mnist
from tensorflow.keras import models
from tensorflow.keras import layers
import matplotlib.pyplot as plt
import numpy as np
import h5py


# load data -------------------------------------------------
data_dir = "/home/kirana/cybertraining2020_team5/research/DATA"

num = 4000 # number of profiles
nchan = 3 # number of channels

# get spatial dimension
fname = data_dir+"/Out1/profile_%05d.hdf5" %(1)
hf = h5py.File(fname, 'r')
spatial = hf.get('x')
x = np.array(spatial)
x_size = np.size(np.array(spatial))
hf.close()

data = np.empty((num, x_size, nchan+1), dtype=float)
for i in range(0, num):
    fname = data_dir+"/Out1/profile_%05d.hdf5" %(i+1)
    hf = h5py.File(fname, 'r')
    data[i, :, 0] = np.array(hf.get("red"))
    data[i, :, 1] = np.array(hf.get("grn"))
    data[i, :, 2] = np.array(hf.get("blu"))
    data[i, :, 3] = np.array(hf.get("tau"))
    hf.close()

np.random.seed(4);
np.random.shuffle(data)

# allocate for training and testing -------------------------
trp = 0.95; # percent to train on

train_size = int(trp * num)
test_size = num - train_size

# spatial slicing dimensions
os = 8 # output slice size
halo = 2; # edge cells
ts = os + halo * 2 # total input slice size including halos
ks = int(ts / 2) # kernel size is half of total size
l2r = 3; # number of slices

print('train_images = ', data[0:5,:,0])
# CNN over domain -------------------------------------------
predictions = np.empty((test_size,os*l2r))
for i in range(0, l2r):

    train_images = data[0:train_size, i*os:i*os+ts, 0:nchan]
    test_images = data[train_size:num, i*os:i*os+ts, 0:nchan]
    train_labels = data[0:train_size, i*os+halo:i*os+halo+os, nchan]
    test_labels = data[train_size:num, i*os+halo:i*os+halo+os, nchan]

    train_shape = train_images.shape
    test_shape = test_images.shape
    print('train_images = ', train_images[0:5,:,:])
    print('test_images = ', test_images[0:5,:,:])
    print('train_labels = ', train_labels[0:5,:])
    print('test_labels = ', test_labels[0:5,:])

    model_m = models.Sequential()
    model_m.add(layers.Conv1D(nchan*100, kernel_size=ks, activation='relu', input_shape=(ts, nchan)))
    model_m.add(layers.Conv1D(4, kernel_size=1))
    model_m.add(layers.Dropout(0.5))
    model_m.add(layers.Flatten())
    model_m.add(layers.Dense(os, activation="linear"))
    #print(model_m.summary())

    # build network ---------------------------------------------

    ltype = 'mean_squared_error'
    bsize = 1024
    eps = 30

    model_m.compile(optimizer='adam',
                    loss=ltype,
                    metrics=[ltype])

    history = model_m.fit(train_images,
                          train_labels,
                          batch_size = bsize,
                          epochs = eps)

    # test and predict ------------------------------------------
    results = model_m.evaluate(x=test_images, y=test_labels)

    predictions[:,i*os:(i+1)*os] = model_m.predict(test_images)


# plotting ------------------------------------------------------
fgnm = "./plots/predict"+"_"+str(ltype)+"_"+str(bsize)+"_"+str(eps)
fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True)

# plots the last three test images
ax = axs[0]
ax.set_title(r"COT")
ax.plot(x[halo:os*l2r+halo],data[num-3,halo:os*l2r+halo,3].T,alpha=0.5,color="blue")
ax.plot(x[halo:os*l2r+halo],predictions[test_size-3,:],alpha=0.5,color="green")
ax.legend(["True", "Predicted"])
ax.set_ylabel(r"case 1")

ax = axs[1]
ax.plot(x[halo:os*l2r+halo],data[num-2,halo:os*l2r+halo,3].T,alpha=0.5,color="blue")
ax.plot(x[halo:os*l2r+halo],predictions[test_size-2,:],alpha=0.5,color="green")
ax.legend(["True", "Predicted"])
ax.set_ylabel(r"case 2")

ax = axs[2]
ax.plot(x[halo:os*l2r+halo],data[num-1,halo:os*l2r+halo,3].T,alpha=0.5,color="blue")
ax.plot(x[halo:os*l2r+halo],predictions[test_size-1,:],alpha=0.5,color="green")
ax.legend(["True", "Predicted"])
ax.set_ylabel(r"case 3")
ax.set_xlabel('X [km]')

plt.savefig(fgnm+".png",dpi=200,bbox_inches='tight')
















