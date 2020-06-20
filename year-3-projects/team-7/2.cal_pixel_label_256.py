'''
This step is to generate 256 by 256 images/data: mask, predictor, figure, composite, and landtype
Input:
    npz files in cal_label folder
    viirs granules in viirs_data folder, after downloading using the download_VIIRS_....sh
    downsized collocated files: e.g., files_north_atlantic.csv
Output:
    This step is to generate 256 by 256 data: mask, predictor, figure, composite, and landtype
    Also generates a records.csv file that records the subset name and how many pixels are dust pixels
'''

import pandas as pd
import numpy as np
from os import listdir, path, mkdir
from pyresample.geometry import AreaDefinition
from satpy import Scene  # need to install satpy library, cannot run on taki
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import from_levels_and_colors

# input folders/data:
project_folder = '/umbc/xfs1/cybertrn/cybertraining2020/team7/research/VIIRS-SIPS/'
cal_label_folder = project_folder+'cal_label/'
viirs_folder = project_folder+'viirs_data/'
input_csv = project_folder+'files_asia_spring.csv'

# output folder:
root = project_folder+'subset_256/asia_spring/'  # for a different spatiotemporal region, new another folder, and change to that folder
predictor_folder = root + 'predictor/'
mask_folder = root + 'mask/'
figure_folder = root + 'figure/'
composite_folder = root + 'composite/'
full_composite = root + 'full_composite/'
lc_folder = root + 'landtype/'

# new a folder if not exists
for folder in [root, predictor_folder, mask_folder, figure_folder, composite_folder, full_composite, lc_folder]:
    if not path.exists(folder):
        mkdir(folder)

matrix_size = 256
stride = 256
bands = ['M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 'M10', 'M11', 'M12', 'M13', 'M14', 'M15', 'M16',
         'satellite_azimuth_angle', 'satellite_zenith_angle', 'solar_azimuth_angle', 'solar_zenith_angle']
cmap, norm = from_levels_and_colors([-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
                                    ['grey', 'white', 'orange', 'green', 'blue', 'purple', 'yellow'])

def fnmatch(fname):
    viirs_files = [file for file in listdir(viirs_folder) if file.endswith(".nc")]
    julian_day = fname[3:10]
    time = fname[-4:]
    VNP02file = [vf for vf in viirs_files if 'VNP02MOD.A' + julian_day + '.' + time in vf][0]
    VNP03file = [vf for vf in viirs_files if 'VNP03MOD.A' + julian_day + '.' + time in vf][0]
    return VNP02file, VNP03file

def checkOverlap(lon, lat, minx, maxx, miny, maxy):
    # check if any of the lon, lat pair is inside the bounding box
    check = 0
    for k in range(lat.shape[0]):
        if lon[k, 0] <= maxx and lon[k, 0] >= minx and lat[k, 0] <= maxy and lat[k, 0] >= miny:
            check = 1
    if check == 0:
        return False
    else:
        return True

df = pd.read_csv(input_csv, index_col=False, header=None, delimiter=' ')
for i in range(20, df.shape[0]):
    fname = df.loc[i, 0].split('.')[-3]
    print(fname)

    cal_label = np.load(cal_label_folder + fname[3:] + '.npz')
    lon = cal_label['lon']
    lat = cal_label['lat']
    CALIOP_Alay_Aerosol_Type_Mode = cal_label['CALIOP_Alay_Aerosol_Type_Mode']
    # N/A: not applicable, 1: clean marine, 2: dust, 3: polluted continental, 4: clean continental, 5: polluted dust, 6:smoke
    CALIOP_N_Clay_5km = cal_label['CALIOP_N_Clay_5km'][:]
    CALIOP_N_Clay_1km = cal_label['CALIOP_N_Clay_1km'][:]
    IGBP_SurfaceType = cal_label['IGBP_SurfaceType'][:]

    if (2 in CALIOP_Alay_Aerosol_Type_Mode):
        minlon, maxlon, minlat, maxlat = df.loc[i, 1:]
        fn1, fn2 = fnmatch(fname)

        # get all available bands
        VNP02 = Dataset(viirs_folder + fn1)
        available_bands = [key for key in VNP02['observation_data'].variables.keys() if len(key) == 3] + [
            'satellite_azimuth_angle', 'satellite_zenith_angle', 'solar_azimuth_angle', 'solar_zenith_angle']

        # load available bands in scene
        scn = Scene(filenames=[viirs_folder + fn1, viirs_folder + fn2],
                    reader='viirs_l1b')  # load VNP02 and VNP03 files together
        scn.load(available_bands + ['dust'])
        dst_area = AreaDefinition('crop_area', 'crop_area', 'crop_latlong', {'proj': 'latlong'},
                                  (maxlon - minlon) / 0.0075, (maxlat - minlat) / 0.0075,
                                  [minlon, minlat, maxlon, maxlat])
        local_scn = scn.resample(dst_area)

        if (local_scn[available_bands[0]].shape[0] > matrix_size) and (
                local_scn[available_bands[0]].shape[1] > matrix_size):
            if not path.exists(full_composite + fname[3:] + '_dust.png'):
                local_scn.save_dataset('dust', full_composite + fname[3:] + '_dust.png')
            composite_image = mpimg.imread(full_composite + fname[3:] + '_dust.png')

            try:
                scn.load(['true_color_raw'])
                local_scn = scn.resample(dst_area)
                local_scn.save_dataset('true_color_raw', full_composite + fname[3:] + '_true_color.png')
                tc_image = mpimg.imread(full_composite + fname[3:] + '_true_color.png')
            except:
                print('no true color')

            category = np.zeros(len(lat))
            for k in range(lat.shape[0]):
                if ((CALIOP_N_Clay_5km[k, 0] == 0) & (CALIOP_N_Clay_1km[k, 0] == 0) & (
                        np.unique(CALIOP_Alay_Aerosol_Type_Mode[:, k]).size == 2) & (
                        2 in CALIOP_Alay_Aerosol_Type_Mode[:, k])):
                    # category 1: dust (dust only no cloud no other aerosols)
                    category[k] = 1
                    continue
                if ((CALIOP_N_Clay_5km[k, 0] == 0) & (CALIOP_N_Clay_1km[k, 0] == 0) & (
                        np.unique(CALIOP_Alay_Aerosol_Type_Mode[:, k]).size == 2) & (
                        5 in CALIOP_Alay_Aerosol_Type_Mode[:, k])):
                    # category 2: pure polluted dust
                    category[k] = 2
                    continue
                if ((CALIOP_N_Clay_5km[k, 0] == 0) & (CALIOP_N_Clay_1km[k, 0] == 0) & (
                        np.unique(CALIOP_Alay_Aerosol_Type_Mode[:, k]).size == 3) & (
                        5 in CALIOP_Alay_Aerosol_Type_Mode[:, k]) & (
                        2 in CALIOP_Alay_Aerosol_Type_Mode[:, k])):
                    # category 3: dust with polluted dust
                    category[k] = 3
                    continue
                if ((CALIOP_N_Clay_5km[k, 0] == 0) & (CALIOP_N_Clay_1km[k, 0] == 0) & (
                        np.unique(CALIOP_Alay_Aerosol_Type_Mode[:, k]).size > 2) & ((
                        2 in CALIOP_Alay_Aerosol_Type_Mode[:, k])) or (
                        5 in CALIOP_Alay_Aerosol_Type_Mode[:, k])):
                    # category 4: dust or polluted dust with other aerosols
                    category[k] = 4
                    continue
                if ((CALIOP_N_Clay_5km[k, 0] == 0) & (CALIOP_N_Clay_1km[k, 0] == 0) & (
                        np.unique(CALIOP_Alay_Aerosol_Type_Mode[:, k]).size > 2) & (
                        2 not in CALIOP_Alay_Aerosol_Type_Mode[:, k]) & (
                        5 not in CALIOP_Alay_Aerosol_Type_Mode[:, k])):
                    # category 5: other aerosol only
                    category[k] = 5
                    continue

            for mi in range(0, local_scn[available_bands[0]].shape[0] - matrix_size, stride):
                for mj in range(0, local_scn[available_bands[0]].shape[1] - matrix_size, stride):
                    x = local_scn[available_bands[0]][mi:mi + matrix_size, mj:mj + matrix_size].x.values
                    y = local_scn[available_bands[0]][mi:mi + matrix_size, mj:mj + matrix_size].y.values

                    if checkOverlap(lon, lat, x.min(), x.max(), y.min(), y.max()):
                        # if any of the lon, lat pair is inside the 128 by 128 image, then write to outputs

                        # save mask/calipso dust into another npy array
                        mask = np.zeros((len(x), len(y)))
                        mask.fill(-1)
                        surface = np.zeros((len(x), len(y)))
                        surface.fill(-1)
                        for xi in range(len(x)):
                            for yj in range(len(y)):
                                for k in range(lat.shape[0]):
                                    if lon[k, 0] <= (x[xi] + 0.0075 / 2) and lon[k, 0] >= (x[xi] - 0.0075 / 2) and lat[
                                        k, 0] <= (y[yj] + 0.0075 / 2) and lat[k, 0] >= (y[yj] - 0.0075 / 2):
                                        mask[xi, yj] = category[k]
                                        surface[xi, yj] = IGBP_SurfaceType[k]
                        print('Dust exist or not:', mask[mask == 1].size)

                        if mask[mask == 1].size > 0:
                            if not path.exists(mask_folder + fname[3:] + '_' + str(mi) + '_' + str(mj) + '.npy'):
                                np.save(mask_folder + fname[3:] + '_' + str(mi) + '_' + str(mj), mask)
                                fig = plt.figure()
                                plt.imshow(np.transpose(mask), cmap=cmap, norm=norm)
                                plt.colorbar(ticks=np.arange(-1, 6))
                                plt.savefig(figure_folder + fname[3:] + '_' + str(mi) + '_' + str(mj))
                                plt.close(fig)

                                with open(root + 'records.csv', 'a') as file:
                                    file.write(fname[3:] + '_' + str(mi) + '_' + str(mj) + ',' + str(
                                        mask[mask == 1].size) + '\n')

                            if not path.exists(lc_folder + fname[3:] + '_' + str(mi) + '_' + str(mj) + '.npy'):
                                np.save(lc_folder + fname[3:] + '_' + str(mi) + '_' + str(mj), surface)

                            # save all available bands into a npy array
                            if not path.exists(predictor_folder + fname[3:] + '_' + str(mi) + '_' + str(mj) + '.npy'):
                                predictors = np.array([])
                                for band in bands:
                                    if band in available_bands:
                                        matrix_data = local_scn[band][mi:mi + matrix_size, mj:mj + matrix_size].values
                                    else:
                                        matrix_data = np.full((matrix_size, matrix_size), np.nan)
                                    if predictors.shape[0] == 0:
                                        predictors = matrix_data
                                    else:
                                        predictors = np.dstack((predictors, matrix_data))
                                print('Predictor size:', predictors.shape)
                                np.save(predictor_folder + fname[3:] + '_' + str(mi) + '_' + str(mj), predictors)

                            # save dust composite image
                            if not path.exists(
                                    composite_folder + fname[3:] + '_' + str(mi) + '_' + str(mj) + '_dust.png'):
                                cropped = composite_image[mi:mi + matrix_size, mj:mj + matrix_size, :]
                                plt.imsave(composite_folder + fname[3:] + '_' + str(mi) + '_' + str(mj) + '_dust.png',
                                            cropped)

                            # save true color composite image
                            if path.exists(full_composite + fname[3:] + '_true_color.png'):
                                if not path.exists(composite_folder + fname[3:] + '_' + str(mi) + '_' + str(
                                        mj) + '_true_color.png'):
                                    cropped_tc = tc_image[mi:mi + matrix_size, mj:mj + matrix_size, :]
                                    plt.imsave(composite_folder + fname[3:] + '_' + str(mi) + '_' + str(
                                        mj) + '_true_color.png',cropped_tc)
