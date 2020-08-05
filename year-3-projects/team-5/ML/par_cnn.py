import tensorflow.keras
from tensorflow.keras.datasets import fashion_mnist
from tensorflow.keras import models
from tensorflow.keras import layers
import matplotlib.pyplot as plt
import numpy as np
import h5py
from mpi4py import MPI
from sklearn.utils import shuffle


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

data_dir = "./../DATA"

channels = ["red", "grn", "blu"]
num = 4000 # number of profiles
nchan = len(channels) # number of channels
trp = 0.95; # percent to train on

os = 8 # output size
halo = 2; # halo around output
ts = os + halo * 2 # input size
ks = int(ts / 2) # input size

nslice = 500; # number of spatial slices to run
l_slice = int(nslice / size)
remainder = nslice % size
ns = np.full((size, 1), l_slice)
for i in range(0, remainder):
    ns[i] = ns[i] + 1
if (rank < remainder):
    l_slice = l_slice + 1

# load data -------------------------------------------------
if rank == 0:

    # get spatial dimension
    fname = data_dir+"/Out1/profile_%05d.hdf5" %(1)
    hf = h5py.File(fname, 'r')
    x = np.array(hf.get('x'))
    hf.close()

    # set up data arrays for host process
    ch1_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    ch2_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    ch3_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    out_data = np.empty((nslice * os, num), dtype='f')

    for i in range(0, num):

        # load data
        fname = data_dir+"/Out1/profile_%05d.hdf5" %(i+1)
        hf = h5py.File(fname, 'r')
        
        ch1_data[:, i] = np.array(hf.get("red"))[:nslice*os+halo*2]
        ch2_data[:, i] = np.array(hf.get("grn"))[:nslice*os+halo*2]
        ch3_data[:, i] = np.array(hf.get("blu"))[:nslice*os+halo*2]
        out_data[:, i] = np.array(hf.get("tau"))[halo:nslice*os+halo]

        hf.close()

    # randomly shuffle (obviously needs a better way to do this)
    np.random.seed(4);
    np.random.shuffle(ch1_data.T)
    np.random.seed(4);
    np.random.shuffle(ch2_data.T)
    np.random.seed(4);
    np.random.shuffle(ch3_data.T)
    np.random.seed(4);
    np.random.shuffle(out_data.T)

else:

    ch1_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    ch2_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    ch3_data = np.empty((nslice * os + halo * 2, num), dtype='f')
    out_data = np.empty((nslice * os, num), dtype='f')

comm.Barrier()

train_size = int(trp * num)
test_size = num - train_size

ch1_recvbuf = np.empty((os*l_slice, num), dtype='f')
ch2_recvbuf = np.empty((os*l_slice, num), dtype='f')
ch3_recvbuf = np.empty((os*l_slice, num), dtype='f')
out_recvbuf = np.empty((os*l_slice, num), dtype='f')

n_elems = ns * os * num

split_disp = np.zeros((size, 1))
split_disp[0] = halo * num
for i in range(1, size):
    split_disp[i] = split_disp[i - 1] + n_elems[i - 1]


comm.Scatterv([ch1_data[0], n_elems, split_disp, MPI.FLOAT], ch1_recvbuf, root = 0)
comm.Scatterv([ch2_data[0], n_elems, split_disp, MPI.FLOAT], ch2_recvbuf, root = 0)
comm.Scatterv([ch3_data[0], n_elems, split_disp, MPI.FLOAT], ch3_recvbuf, root = 0)
split_disp_o = split_disp - halo * num
comm.Scatterv([out_data[0], n_elems, split_disp_o, MPI.FLOAT], out_recvbuf, root = 0)

data_images = np.empty((num, os * l_slice + 2 * halo, nchan), dtype='f')
data_images[:,halo:l_slice*os+halo,0] = ch1_recvbuf[:,:num].T
data_images[:,halo:l_slice*os+halo,1] = ch2_recvbuf[:,:num].T
data_images[:,halo:l_slice*os+halo,2] = ch3_recvbuf[:,:num].T
data_labels = out_recvbuf[:,:].T

if rank == 0:
    beg1_local = np.empty((halo, num), dtype='f')
    beg2_local = np.empty((halo, num), dtype='f')
    beg3_local = np.empty((halo, num), dtype='f')
    eeg1_local = np.empty((halo, num), dtype='f')
    eeg2_local = np.empty((halo, num), dtype='f')
    eeg3_local = np.empty((halo, num), dtype='f')
else:
    beg1_local = np.empty((halo, num), dtype='f')
    beg2_local = np.empty((halo, num), dtype='f')
    beg3_local = np.empty((halo, num), dtype='f')
    eeg1_local = np.empty((halo, num), dtype='f')
    eeg2_local = np.empty((halo, num), dtype='f')
    eeg3_local = np.empty((halo, num), dtype='f')

split_disp_b = split_disp - halo * num
comm.Scatterv([ch1_data[0], halo*num, split_disp_b, MPI.FLOAT], beg1_local, root = 0)
comm.Scatterv([ch2_data[0], halo*num, split_disp_b, MPI.FLOAT], beg2_local, root = 0)
comm.Scatterv([ch3_data[0], halo*num, split_disp_b, MPI.FLOAT], beg3_local, root = 0)

split_disp_e = split_disp + n_elems
comm.Scatterv([ch1_data[0], halo*num, split_disp_e, MPI.FLOAT], eeg1_local, root = 0)
comm.Scatterv([ch2_data[0], halo*num, split_disp_e, MPI.FLOAT], eeg2_local, root = 0)
comm.Scatterv([ch3_data[0], halo*num, split_disp_e, MPI.FLOAT], eeg3_local, root = 0)

data_images[:,0:halo,0] = beg1_local[:,:num].T
data_images[:,0:halo,1] = beg2_local[:,:num].T
data_images[:,0:halo,2] = beg3_local[:,:num].T
data_images[:,-halo:,0] = eeg1_local[:,:num].T
data_images[:,-halo:,1] = eeg2_local[:,:num].T
data_images[:,-halo:,2] = eeg3_local[:,:num].T


# run ML model -----------------------------------------------------------
l_predict = np.empty((os*l_slice, test_size), dtype='f')
for i in range(0, l_slice):

    train_images = data_images[0:train_size, i*os:i*os+ts, 0:nchan]
    test_images = data_images[train_size:num, i*os:i*os+ts, 0:nchan]
    train_labels = data_labels[0:train_size, i*os:i*os+os]
    test_labels = data_labels[train_size:num, i*os:i*os+os]

    train_shape = train_images.shape
    test_shape = test_images.shape

    model_m = models.Sequential()
    model_m.add(layers.Conv1D(nchan*100, kernel_size=ks, activation='relu', input_shape=(ts, nchan)))
    model_m.add(layers.Conv1D(4, kernel_size=1))
    model_m.add(layers.Dropout(0.5))
    model_m.add(layers.Flatten())
    model_m.add(layers.Dense(os, activation="linear"))

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
    l_predict[i*os:(i+1)*os,:] = model_m.predict(test_images).T
comm.Barrier()

# move data back to proc 0 --------------------------------------
predictions = np.empty((os*nslice, test_size), dtype='f')
n_elems = ns * os * test_size
split_disp = np.zeros((size, 1))
split_disp[0] = 0
for i in range(1, size):
    split_disp[i] = split_disp[i - 1] + n_elems[i - 1]
comm.Gatherv(l_predict, [predictions, nslice*os*test_size, split_disp, MPI.FLOAT], root=0)
comm.Barrier()

# plotting ------------------------------------------------------
if rank == 0:

    fgnm = "./plots/par_predict"+"_"+str(ltype)+"_"+str(bsize)+"_"+str(eps)
    fig, axs = plt.subplots(nrows=3,ncols=1,sharex=True)

    # plots the last three test images
    ax = axs[0]
    ax.plot(x[halo:os*nslice+halo],out_data[:,num-3],color="blue")
    ax.plot(x[halo:os*nslice+halo],predictions[:,test_size-3].T,alpha=0.5,color="green")
    ax.legend(["True", "Predicted"])
    ax.set_ylabel(r"case 1")

    ax = axs[1]
    ax.plot(x[halo:os*nslice+halo],out_data[:,num-2],color="blue")
    ax.plot(x[halo:os*nslice+halo],predictions[:,test_size-2].T,alpha=0.5,color="green")
    ax.legend(["True", "Predicted"])
    ax.set_ylabel(r"case 2")

    ax = axs[2]
    ax.plot(x[halo:os*nslice+halo],out_data[:,num-1],color="blue")
    ax.plot(x[halo:os*nslice+halo],predictions[:,test_size-1].T,alpha=0.5,color="green")
    ax.legend(["True", "Predicted"])
    ax.set_ylabel(r"case 3")
    ax.set_xlabel('X [km]')

    plt.savefig(fgnm+".png",dpi=200,bbox_inches='tight')


