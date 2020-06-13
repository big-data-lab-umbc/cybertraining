'''
Use three different clustering methods to cluster all images within the test set.
Input:
    change the model_index number for a different clustering method
    test set folder, e.g., viirs_north_atlantic
    K number: how many clusters
Output:
    csv files that record the accuracy metrics along the track of calipso, and over different surface types
    generate box plots of accuracy metrics using different clustering methods and over different surface types
'''

import numpy as np
from os import listdir
from sklearn.cluster import KMeans #K-Means Clustering
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import pandas as pd
import matplotlib.pyplot as plt
from sklearn_extra.cluster import KMedoids
from fcmeans import FCM

################################################ Step 1: Load the data ################################################

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
accuracy_record = root+'accuracy_k' + str(4) + '_'+selected_models[model_index][1]+'.csv'

land_types = {-1: 'N/A', 1:'Evergreen Needleleaf Forests', 2:'Evergreen Broadleaf Forests',3:'Deciduous Needleleaf Forests',
              4: 'Deciduous Broadleaf Forests', 5: 'Mixed Forests', 6:'Closed Shrublands', 7:'Open Shrublands', 8:'Woody Savannas',
              9:'Savannas', 10:'Grasslands', 11:'Permanent Wetlands', 12:'Croplands', 13:' Urban and Built-Up Lands',
              14:'Cropland/Natural Vegetation Mosaics', 15:'Snow and Ice', 16:'Barren', 17:'Water Bodies'}
land_type_values = list(land_types.keys())
land_type_values.sort()

with open(accuracy_record, 'a') as file:
    file.write('fname,accuracy,precision,recall,f1-score,ocean_accuracy,ocean_precision,ocean_recall,ocean_f1-score,'
               'barren_accuracy,barren_precision,barren_recall,barren_f1-score,other_accuracy,other_precision,other_recall,other_f1-score\n')

fnames = [file.split('.')[0] for file in listdir(mask_folder)]
for f in range(len(fnames)):
    # load viirs all bands
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

        ############################################## Step 2: clustering ##############################################
        model = selected_models[model_index][0]
        # Fit our model
        model.fit(vectorized)
        # Find which cluster each data-point belongs to
        labels = model.predict(vectorized)
        label_image = labels.reshape(predictor.shape[0:2]) #Reshape the clusters back to image size
        # print(label_image[mask == 1])
        counts = np.bincount(label_image[mask == 1])
        try:
            print(np.argmax(counts)) # find the majority cluster along the track

            # Assuming that the majority cluster on the mask track is the dust class
            dust_label = np.argmax(counts)
            dust = label_image.copy()
            dust[dust != dust_label] = 0

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

            with open(accuracy_record, 'a') as file:
                file.write(fnames[f]+',' + str(accuracy_score(y_true, y_pred)) + ',' + str(precision_score(y_true, y_pred))+','
                           + str(recall_score(y_true, y_pred))+',' + str(f1_score(y_true, y_pred))+ ','
                           + str(a1)+ ',' + str(a2) + ',' + str(a3)+ ',' + str(a4) + ','
                           + str(b1)+ ',' + str(b2) + ',' + str(b3)+ ',' + str(b4) + ','
                           + str(c1)+ ',' + str(c2) + ',' + str(c3)+ ',' + str(c4) +'\n')
        except ValueError:
            pass
    except FileNotFoundError:
        pass

df = pd.read_csv(accuracy_record, index_col=False)
fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df['accuracy'], df['precision'], df['recall'], df['f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of '+selected_models[model_index][1]+' accuracy metrics with 1037 images')
plt.show()

# test accuracy over the ocean/water body
df_ocean = df[df['ocean_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_ocean['ocean_accuracy'], df_ocean['ocean_precision'], df_ocean['ocean_recall'], df_ocean['ocean_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of '+selected_models[model_index][1]+' over waterbodies with 249 images')
plt.show()

# test accuracy over barren
df_barren = df[df['barren_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_barren['barren_accuracy'], df_barren['barren_precision'], df_barren['barren_recall'], df_barren['barren_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of '+selected_models[model_index][1]+' over barren with 648 images')
plt.show()

# test accuracy over other surface types
df_other = df[df['other_accuracy']!=-1]

fig = plt.figure()
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot([df_other['other_accuracy'], df_other['other_precision'], df_other['other_recall'], df_other['other_f1-score']])
ax.set_xticklabels(['Accuracy', 'Precision', 'Recall', 'F1-score'])
plt.title('Boxplot of '+selected_models[model_index][1]+' over other surface types with 207 images')
plt.show()