'''
Use K-means clustering to cluster one image within the test set.
Input:
    test set folder, e.g., viirs_north_atlantic
    K number: how many clusters
Output:
    for this test image: K-means clusters, ground truth/calipso track dust categories, landtype plot
'''

import math
import pandas as pd
from os import listdir
import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, pdist
from sklearn.cluster import KMeans
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix
from matplotlib.colors import from_levels_and_colors

root = 'I:\\viirs_north_africa_summer\\'
predictor_folder = root + 'predictor\\'
mask_folder = root + 'mask\\'
figure_folder = root + 'figure\\'
composite_folder = root + 'composite\\'
full_composite = root + 'full_composite\\'
lc_folder = root + 'landtype\\'

land_types = {-1: 'N/A', 0:'N/A', 1:'Evergreen Needleleaf Forests', 2:'Evergreen Broadleaf Forests',3:'Deciduous Needleleaf Forests',
              4: 'Deciduous Broadleaf Forests', 5: 'Mixed Forests', 6:'Closed Shrublands', 7:'Open Shrublands', 8:'Woody Savannas',
              9:'Savannas', 10:'Grasslands', 11:'Permanent Wetlands', 12:'Croplands', 13:' Urban and Built-Up Lands',
              14:'Cropland/Natural Vegetation Mosaics', 15:'Snow and Ice', 16:'Barren', 17:'Water Bodies'}

figure_size = 5

def plot_segmentation(labels_image, title):
    plt.figure(figsize=(figure_size, figure_size))
    cmap = plt.get_cmap('RdBu', np.max(labels_image) - np.min(labels_image) + 1)
    im = plt.imshow(labels_image, cmap=cmap, vmin=np.min(labels_image) - .5, vmax=np.max(labels_image) + .5)
    cax = plt.colorbar(im, ticks=np.arange(np.min(labels_image), np.max(labels_image) + 1), fraction=0.05, pad=0.05)
    plt.title(title), plt.xticks([]), plt.yticks([])
    plt.show()

def plot_img(img, title):
    plt.figure(figsize=(figure_size,figure_size))
    plt.imshow(img)
    plt.title(title), plt.xticks([]), plt.yticks([])
    plt.show()

def plot_dust(plot_label, title):
    cmap, norm = from_levels_and_colors(np.linspace(-1.5, 5.5, num=8),['grey','white','orange','green','blue','purple','yellow'])
    plt.figure()
    plt.imshow(plot_label, cmap=cmap, norm=norm)
    plt.colorbar(ticks=np.arange(-1, 6))
    plt.title(title), plt.xticks([]), plt.yticks([])
    plt.show()

def plot_lc(landtype):
    plt.figure()
    min, max = np.min(list(land_types.keys())), np.max(list(land_types.keys()))
    cmap, norm = from_levels_and_colors(np.linspace(-1.5, 17.5, num=19),
                                        ['grey','grey','green','green','green','green','green','green','green','green','green','green','green','brown', 'green', 'white', 'yellow','blue'])
    im = plt.imshow(landtype, cmap=cmap, vmin=min - .5, vmax=max + .5)
    ticks = list(land_types.keys())
    ticks.sort()
    cbar = plt.colorbar(im, ticks=ticks, fraction=0.05, pad=0.05)
    cbar.ax.set_yticklabels([land_types[i] for i in ticks])
    plt.title('landtype'), plt.xticks([]), plt.yticks([])
    plt.show()

def dust_predictors():
    # df = pd.read_csv(root + 'records.csv', header=None, index_col=False)
    # df = df[df[1] >= 100]
    # fnames = df[0].values
    # fnames = [file+'.npy' for file in fnames]

    fnames = [file for file in listdir(predictor_folder)]

    aggr_predictor = []
    aggr_mask = []
    # load viirs all bands
    for fname in fnames:
        predictor_file = predictor_folder+fname
        predictor = np.load(predictor_file)
        predictor[np.isnan(predictor)] = 0
        # print(predictor.shape)
        predictor = predictor[:,:,:16]

        vectorized = predictor.reshape((-1, predictor.shape[2]))
        vectorized = np.float32(vectorized)
        aggr_predictor.append(vectorized)

        # load calipso dust mask
        mask_file = mask_folder + fname
        mask = np.transpose(np.load(mask_file))
        mask_vectorized = mask.reshape((-1))
        aggr_mask.append(mask_vectorized)

    aggr_predictor = np.concatenate(aggr_predictor, axis = 0)
    aggr_mask = np.concatenate(aggr_mask, axis=0)

    return aggr_predictor[aggr_mask == 1]

# African dust
# fname = '2014152t0936_256_0'     # good with pure dust (pure dust only)
# fname = '2014152t1112_1024_256'  # good with categories 1 and 3
fname = '2014152t1254_1536_256'
#fname = '2014152t1436_1024_256'
#fname = '2014152t1436_1280_256'
#fname = '2014152t1436_1536_256'

# North Atlantic
# fname = '2014187t0548_2304_0'
# fname = '2014360t0324_1792_0'
# fname = '2014245t1536_1024_256'
# fname = '2014211t0500_2048_0'
# fname = '2014229t1718_1536_256'
# fname = '2014275t0500_2304_0'
# fname = '2014240t1530_768_0'
# fname = '2014240t1530_1024_256'
# fname = '2014234t1724_512_0'

# Asian dust
# fname = '2014136t0430_256_0'
    #'2014136t0430_256_0'
    #'2014136t0430_1536_256'
    #'2014131t0606_1280_256'
    #'2014072t1830_1280_256'
    #'2014072t1830_1536_256'
    #'2014072t1830_1792_0'
    #'2014147t0606_256_0'

# load viirs all bands
predictor_file = predictor_folder+fname+'.npy'
predictor = np.load(predictor_file)
predictor[np.isnan(predictor)] = 0
# print(predictor.shape)
predictor = predictor[:,:,:16]

# load calipso dust mask
mask_file = mask_folder+fname+'.npy'
mask = np.transpose(np.load(mask_file))
# print(mask.shape)

# Surface type
land_file = lc_folder + fname + '.npy'
land = np.transpose(np.load(land_file))
print([land_types[type] for type in np.unique(land)])

# load dust composite image
img_file = composite_folder+fname+'_dust.png'
print(img_file)
original_image = cv2.imread(img_file)
img = cv2.cvtColor(original_image,cv2.COLOR_BGR2RGB)

vectorized = predictor.reshape((-1,predictor.shape[2]))
vectorized = np.float32(vectorized)

#Initialize our model
K = 4
model = KMeans(n_clusters=K)
#Fit our model
model.fit(vectorized)
#Find which cluster each data-point belongs to
labels = model.predict(vectorized)
#Reshape the clusters back to image size
labels_image = labels.reshape(predictor.shape[0:2])
centroids = model.cluster_centers_
print(labels_image.shape)

plot_img(img,'Dust composite')
plot_segmentation(labels_image, 'Segmentation with K='+str(K))

# Assuming that the majority cluster on the mask track is the pure dust class
# counts = np.bincount(plot_label[mask == 1])
# dust_label = np.argmax(counts)
# # print(dust_label)

# improve cluster determination based on dust predictors
dust_profiles = dust_predictors()
mean_distances = []
for i in range(K):
    # Use all pixels to calculate similarity
    mean_distances.append([i, cdist(dust_profiles, predictor[labels_image == i], 'euclidean').mean()])
    # Use centroids to calculate similarity：
    # mean_distances.append([i, np.linalg.norm(dust_profiles.mean(axis=0)-centroids[i])])
mean_distances = np.array(mean_distances)
min_mean_distance = min(mean_distances[:,1])
dust_label = np.where(mean_distances[:,1] == min(mean_distances[:,1]))[0][0]
# check if there are other clusters that have strong similarities with the dust predictors
potentials = []
for i in range(K):
    if i != dust_label:
        if (mean_distances[i, 1] - min_mean_distance) < 10:
            potentials.append(i)
if len(potentials) > 0:
    for p in potentials:
        labels_image[labels_image == p] = dust_label # assign this potential cluster to the original dust cluster

# plot of clusters along track
plot_label = labels_image.copy().astype('Int32')
plot_label[mask < 0] = -1
plot_dust(plot_label, "Kmeans clusters along track")
plot_dust(mask, "Dust mask along track")
plot_lc(land)

y_pred = plot_label[mask >= 0]
y_pred[y_pred != dust_label] = -100
y_pred[y_pred == dust_label] = 100

y_pred[y_pred == -100] = 0
y_pred[y_pred == 100] = 1
# print('y_pred: ', y_pred)

y_true = mask[mask >= 0].flatten()
y_true[y_true == 1] = 1
y_true[y_true != 1] = 0
# print('y_true: ', y_true)

print('Pure dust detection report')
print(classification_report(y_true, y_pred))

print('Pure dust detection accuracy')
print(accuracy_score(y_true, y_pred))

print('Confusion matrix')
print(confusion_matrix(y_true, y_pred))

print('final dust plot')
dust = labels_image.copy()
dust[dust != dust_label] = -100
dust[dust == dust_label] = 100

dust[dust == -100] = 0
dust[dust == 100] = 1

plot_dust(dust, "Dust extent")
