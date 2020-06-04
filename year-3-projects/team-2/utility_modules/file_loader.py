# Third party imports.
import numpy as np
from pyhdf.SD import SD, SDC

# Local application imports.
from vertical_feature_mask import extract_features

# /umbc/xfs1/zzbatmos/common/Data/CALIPSO/CAL_LID_L2_05kmCLay/2007/CAL_LID_L2_05kmCLay-Standard-V4-20.2007-01-11T05-56-51ZN.hdf


def _read_file(filename):
    
    # TODO: Add logging for filename and exception measured.
    
    try:
        
        file = SD(filename, SDC.READ)
        
    except Exception:
        
        file = None
        
        print(filename)
        
    return file


def _generate_files(filenames):
    
    """
    
    A generator function which accepts a pandas.DataFrame of HDF4 filenames
    and converts each filename rowwise into a pyhdf.SD.SD object.
    
    """
    
    for index, row in filenames.iterrows():
        
        files = row.apply(_read_file)
        
        yield files


def _get_data(filenames, datasets):
    
    """
    
    A wrapper generator for the "_generate_files" generator which reads in the
    desired datasets from the loaded files.
    
    """
    
    for files in _generate_files(filenames):
        
        try:
        
            data = {}

            for key, flags in datasets.items():

                for flag in flags:

                    array = np.asarray(files[key].select(flag).get())

                    data[flag] = array

            yield data
            
        except Exception:
            
            yield None
        

def _get_collocated_data(filenames, datasets):
    
    """
    
    A wrapper generator for the "_get_data" generator which selects only those
    values which are collocated.
    
    """
    
    for data in _get_data(filenames, datasets):
        
        # NOTE: Some CALIPSO data flags, such as the Opacity_Flag, exist in the
        #       5km data products but not the 1km data products. By definition,
        #       the 5km product has 1/5th as many data points as does the 1km
        #       product. Consequently, the collocation masking must take this
        #       into consideration. We choose to select only those 5km CALIPSO
        #       columns for which all 1km products are collocated with the
        #       MODIS data.
        if data is not None:
        
            mask_shape_1km = data["Collocation_Flag"].shape
            mask_shape_5km = (data["Collocation_Flag"].size // 5, 5)

            collocation_mask_5km = data["Collocation_Flag"].reshape(mask_shape_5km)
            collocation_mask_5km = collocation_mask_5km.all(axis = 1)
            collocation_mask_1km = np.repeat(collocation_mask_5km, 5)

            for key, value in data.items():

                if data[key].shape[0] == mask_shape_1km[0]:

                    data[key] = value[collocation_mask_1km]

                elif data[key].shape[0] == mask_shape_5km[0]:

                    data[key] = value[collocation_mask_5km]

            yield data
            
        else:
            
            continue
        
        
def _preprocess_data(filenames, datasets):
    
    """
    
    A wrapper generator for the "_get_collocated_data" generator which
    performs some basic preprocessing on the collocated data, such as
    extracting information from the feature classification flags.
    
    """
    
    # TODO: Check that keys exist first.
    
    for index, data in enumerate(_get_collocated_data(filenames, datasets)):
        
        if data is not None:
            
            files = filenames.iloc[index]
            files = {**files.to_dict(), **{"timestamp" : files.name}}
            
            data = {**data, **files}

            data["Feature_Classification_Flags"] = \
                extract_features(data["Feature_Classification_Flags"])

            data["Opacity_Flag"] = data["Opacity_Flag"].astype(np.float16)
            data["Opacity_Flag"][data["Opacity_Flag"] == 99.0] = np.nan

            data["MYD06_Cloud_Optical_Thickness"] = data["MYD06_Cloud_Optical_Thickness"].T[0]
            data["MYD06_Cloud_Optical_Thickness"][data["MYD06_Cloud_Optical_Thickness"] == -9999.99] = np.nan
            
            try:
                
                data["MYD06_Cloud_Top_Height_1km"][data["MYD06_Cloud_Top_Height_1km"] == -9999.99] = np.nan
                data["MYD06_Cloud_Top_Height_1km"] = data["MYD06_Cloud_Top_Height_1km"].T[0]
                
            except Exception as e:
                
                print(e)
            
            #data["Column_Optical_Depth_Cloud_532"] = data["Column_Optical_Depth_Cloud_532"].T[0]

            yield data
            
        else:
            
            continue
        
