'''
Use three different clustering methods to cluster one image within the test set.
Input:
    change the model_index number for a different clustering method
    test set folder, e.g., viirs_north_atlantic
    K number: how many clusters
Output:
    for this test image: K-means clusters, ground truth/calipso track dust categories, landtype plot
'''

from sklearn_extra.cluster import KMedoids
from fcmeans import FCM
import math
import pandas as pd
import numpy as np
import cv2 # pip install opencv-python
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import confusion_matrix
from matplotlib.colors import from_levels_and_colors

root = 'I:\\viirs_north_atlantic\\'
predictor_folder = root + 'predictor\\'
mask_folder = root + 'mask\\'
figure_folder = root + 'figure\\'
composite_folder = root + 'composite\\'
full_composite = root + 'full_composite\\'
lc_folder = root + 'landtype\\'

K = 4
selected_models = [
    (KMedoids(metric="euclidean", n_clusters=K), "KMedoids"),
    (KMeans(n_clusters=K), "KMeans"),
    (FCM(n_clusters=K), "Fuzzy C-Means")
]
model_index = 2

# fname = '2014152t1112_0_0'
# fname = '2014152t0936_256_0'     # good with pure dust (pure dust only)
# fname = '2014152t0936_1024_256'  # good with categories 1 and 4
# fname = '2014152t1112_1024_256'  # good with categories 1 and 3
# fname = '2014179t1112_1536_512' # over waterbody


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
    plt.title(title)
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
model = selected_models[model_index][0]
#Fit our model
model.fit(vectorized)
#Find which cluster each data-point belongs to
labels = model.predict(vectorized)
# centroids = model.cluster_centers_
#Reshape the clusters back to image size
labels_image = labels.reshape(predictor.shape[0:2])
print(labels_image.shape)

plot_img(img,'Dust composite')
plot_segmentation(labels_image, 'Segmentation using '+selected_models[model_index][1]+' with K='+str(K))

# plot of clusters along track
plot_label = labels_image.copy().astype('Int32')
plot_label[mask < 0] = -1
plot_dust(plot_label, selected_models[model_index][1]+" clusters along track")
plot_dust(mask, "Dust mask along track")
plot_lc(land)

# Assuming that the majority cluster on the mask track is the pure dust class
counts = np.bincount(plot_label[mask == 1])
dust_label = np.argmax(counts)
# print(dust_label)

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

print(selected_models[model_index][1]+': Pure dust detection report')
print(classification_report(y_true, y_pred))

print(selected_models[model_index][1]+': Pure dust detection accuracy')
print(accuracy_score(y_true, y_pred))

print(selected_models[model_index][1]+': Confusion matrix')
print(confusion_matrix(y_true, y_pred))

print('Final dust plot')
dust = labels_image.copy()
dust[dust != dust_label] = -100
dust[dust == dust_label] = 100

dust[dust == -100] = 0
dust[dust == 100] = 1

plot_dust(dust, "Dust extent based on "+selected_models[model_index][1])