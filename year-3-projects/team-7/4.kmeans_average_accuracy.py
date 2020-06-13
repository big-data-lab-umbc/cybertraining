'''
Use kmeans clustering method to cluster all images within the test set.
Input:
    test set folder, e.g., viirs_north_atlantic
    K number: how many clusters
Output:
    root+'accuracy_k4.csv' csv file that records the accuracy metrics along the track of calipso, and over different surface types
    generate box plots of accuracy metrics using different clustering methods and over different surface types
'''

import numpy as np
from os import listdir
from sklearn.cluster import KMeans #K-Means Clustering
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

################################################ Step 1: Load the data ################################################

root = 'I:\\viirs_north_africa_summer\\'  # for a different spatiotemporal region, new another folder, and change to that folder
predictor_folder = root + 'predictor\\'
mask_folder = root + 'mask\\'
figure_folder = root + 'figure\\'
composite_folder = root + 'composite\\'
full_composite = root + 'full_composite\\'
lc_folder = root + 'landtype\\'

land_types = {-1: 'N/A', 1:'Evergreen Needleleaf Forests', 2:'Evergreen Broadleaf Forests',3:'Deciduous Needleleaf Forests',
              4: 'Deciduous Broadleaf Forests', 5: 'Mixed Forests', 6:'Closed Shrublands', 7:'Open Shrublands', 8:'Woody Savannas',
              9:'Savannas', 10:'Grasslands', 11:'Permanent Wetlands', 12:'Croplands', 13:' Urban and Built-Up Lands',
              14:'Cropland/Natural Vegetation Mosaics', 15:'Snow and Ice', 16:'Barren', 17:'Water Bodies'}
land_type_values = list(land_types.keys())
land_type_values.sort()

K = 10
with open(root+'accuracy_k' + str(K) + '.csv', 'a') as file:
    file.write('fname,accuracy,precision,recall,f1-score,ocean_accuracy,ocean_precision,ocean_recall,ocean_f1-score,'
               'barren_accuracy,barren_precision,barren_recall,barren_f1-score,other_accuracy,other_precision,other_recall,other_f1-score \n')

def dust_predictors():
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

dust_profiles = dust_predictors()

# fnames = [file.split('.')[0] for file in listdir(mask_folder)]
df = pd.read_csv(root+'records.csv', header=None, index_col=False)
df = df[df[1] >= 100]
fnames = df[0].values
for f in range(len(fnames)):    # load viirs all bands
    try:
        predictor_file = predictor_folder+fnames[f]+'.npy'
        predictor = np.load(predictor_file)
        predictor[np.isnan(predictor)] = 0 # fill nan values as 0
        predictor = predictor[:,:,:16] # use only the first 16 bands, no angles used here
        # print(predictor.shape)

        # load calipso dust mask
        mask_file = mask_folder+fnames[f]+'.npy'
        mask = np.transpose(np.load(mask_file))
        # print(mask.shape)

        # Surface type
        land_file = lc_folder + fnames[f]+'.npy'
        land = np.transpose(np.load(land_file))
        try:
            print('Land types', [land_types[type] for type in np.unique(land)])
        except KeyError:
            print(fnames, 'contains strange land type')
            continue
        on_track_land_type = land[mask >= 0]

        # reshape the predictors
        vectorized = predictor.reshape((-1,predictor.shape[2]))
        vectorized = np.float32(vectorized)

        ############################################## Step 2: K-means clustering ##############################################

        kmeans = KMeans(n_clusters=K)
        kmeans.fit(vectorized) #Fit our model
        clusters = kmeans.predict(vectorized) #Find which cluster each data-point belongs to
        label_image = clusters.reshape(predictor.shape[0:2]) #Reshape the clusters back to image size
        # print(label_image[mask == 1])
        counts = np.bincount(label_image[mask == 1])
        try:
            # print(np.argmax(counts)) # find the majority cluster along the track
            #
            # # Assuming that the majority cluster on the mask track is the dust class
            # dust_label = np.argmax(counts)
            # dust = label_image.copy()
            # dust[dust != dust_label] = 0

            # improve cluster determination based on dust predictors

            mean_distances = []
            for i in range(K):
                mean_distances.append([i, cdist(dust_profiles, predictor[label_image == i], 'euclidean').mean()])
            mean_distances = np.array(mean_distances)
            min_mean_distance = min(mean_distances[:, 1])
            dust_label = np.where(mean_distances[:, 1] == min(mean_distances[:, 1]))[0][0]
            # check if there are other clusters that have strong similarities with the dust predictors
            potentials = []
            for i in range(K):
                if i != dust_label:
                    if (mean_distances[i, 1] - min_mean_distance) < 10:
                        potentials.append(i)
            if len(potentials) > 0:
                for p in potentials:
                    label_image[
                        label_image == p] = dust_label  # assign this potential cluster to the original dust cluster

            # Evaluate on track classification accuracy
            y_pred = label_image[mask >= 0]
            y_pred[y_pred != dust_label] = -100
            y_pred[y_pred == dust_label] = 100

            y_pred[y_pred == -100] = 0
            y_pred[y_pred == 100] = 1

            # print(y_pred)
            y_true = mask[mask >= 0].flatten()
            y_true[y_true == 1] = 1
            y_true[y_true != 1] = 0

            print(accuracy_score(y_true, y_pred), precision_score(y_true, y_pred), recall_score(y_true, y_pred), f1_score(y_true, y_pred))

            # ocean accuracy
            if 17 in on_track_land_type:
                ocean_y_pred = y_pred[on_track_land_type==17]
                ocean_y_true = y_true[on_track_land_type==17]
                a1, a2, a3, a4 = accuracy_score(ocean_y_true, ocean_y_pred), precision_score(ocean_y_true, ocean_y_pred), recall_score(ocean_y_true, ocean_y_pred), f1_score(ocean_y_true, ocean_y_pred)
                print(accuracy_score(ocean_y_true, ocean_y_pred), precision_score(ocean_y_true, ocean_y_pred), recall_score(ocean_y_true, ocean_y_pred),
                      f1_score(ocean_y_true, ocean_y_pred))
            else:
                a1, a2, a3, a4 = -1, -1, -1, -1

            # bare land accuracy
            if 16 in on_track_land_type:
                barren_y_pred = y_pred[on_track_land_type==16]
                barren_y_true = y_true[on_track_land_type==16]
                b1, b2, b3, b4 = accuracy_score(barren_y_true, barren_y_pred), precision_score(barren_y_true, barren_y_pred), recall_score(barren_y_true, barren_y_pred), f1_score(barren_y_true, barren_y_pred)
                print(accuracy_score(barren_y_true, barren_y_pred), precision_score(barren_y_true, barren_y_pred), recall_score(barren_y_true, barren_y_pred),
                      f1_score(barren_y_true, barren_y_pred))
            else:
                b1, b2, b3, b4 = -1, -1, -1, -1

            # other:
            if (not 16 in on_track_land_type) and (not 17 in on_track_land_type):
                other_y_pred = y_pred[(on_track_land_type!=16) & (on_track_land_type!=17)]
                other_y_true = y_true[(on_track_land_type != 16) & (on_track_land_type != 17)]
                c1, c2, c3, c4 = accuracy_score(other_y_true, other_y_pred), precision_score(other_y_true, other_y_pred), recall_score(other_y_true, other_y_pred), f1_score(other_y_true, other_y_pred)
                print(accuracy_score(other_y_true, other_y_pred), precision_score(other_y_true, other_y_pred), recall_score(other_y_true, other_y_pred),
                      f1_score(other_y_true, other_y_pred))
            else:
                c1, c2, c3, c4 = -1, -1, -1, -1

            with open(root+'accuracy_k'+str(K)+'.csv', 'a') as file:
                file.write(fnames[f]+',' + str(accuracy_score(y_true, y_pred)) + ',' + str(precision_score(y_true, y_pred))+','
                           + str(recall_score(y_true, y_pred))+',' + str(f1_score(y_true, y_pred))
                           + ','
                           + str(a1)+ ',' + str(a2) + ',' + str(a3)+ ',' + str(a4) + ','
                           + str(b1)+ ',' + str(b2) + ',' + str(b3)+ ',' + str(b4) + ','
                           + str(c1)+ ',' + str(c2) + ',' + str(c3)+ ',' + str(c4)
                           +'\n')
        except ValueError:
            pass
    except FileNotFoundError:
        pass

df = pd.read_csv(root+'accuracy_k'+str(K)+'.csv', index_col=False)
fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df['accuracy'], df['precision'], df['recall'], df['f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of accuracy metrics with '+str(df.shape[0])+' images')
plt.show()

# test accuracy over the ocean/water body
df_ocean = df[df['ocean_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_ocean['ocean_accuracy'], df_ocean['ocean_precision'], df_ocean['ocean_recall'], df_ocean['ocean_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of accuracy metrics over waterbodies with '+str(df_ocean.shape[0])+' images')
plt.show()

# test accuracy over barren
df_barren = df[df['barren_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_barren['barren_accuracy'], df_barren['barren_precision'], df_barren['barren_recall'], df_barren['barren_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of accuracy metrics over barren with '+str(df_barren.shape[0])+' images')
plt.show()

# test accuracy over other surface types
df_other = df[df['other_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_other['other_accuracy'], df_other['other_precision'], df_other['other_recall'], df_other['other_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of accuracy metrics over other surface types with '+str(df_other.shape[0])+' images')
plt.show()