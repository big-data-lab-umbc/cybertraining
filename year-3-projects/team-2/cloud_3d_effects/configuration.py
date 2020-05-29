# Standard library imports.
import os


def get_taki_directories():

    directories = {}

    directories["group_data"] = \
        os.path.abspath("/umbc/xfs1/zzbatmos/common/Data/")

    directories["collocation_files"] = \
        os.path.join(directories["group_data"],
                     "CALIPSO-MODIS-CloudSat")

    directories["collocation_data"] = \
        os.path.join(directories["collocation_files"],
                     "Data")

    directories["collocation_indices"] = \
        os.path.join(directories["collocation_files"],
                     "Index")

    directories["CALIPSO_01km_data"] = \
        os.path.join(directories["group_data"],
                     "CALIPSO/CAL_LID_L2_01kmCLay")
        
#    directories["CALIPSO_05km_data"] = \
#        os.path.join(directories["group_data"],
#                     "CALIPSO/CAL_LID_L2_05kmCLay")
        
    del directories["group_data"]
    del directories["collocation_files"]
        
    return directories


def get_local_directories():
    
    directories = {}

    directories["group_data"] = \
        os.path.abspath("C:/Users/Erick Shepherd/Documents/Data")

    directories["collocation_files"] = \
        os.path.join(directories["group_data"])

    directories["collocation_data"] = \
        os.path.join(directories["collocation_files"],
                     "2007_Data_CALIPSO-MODIS-CloudSat")

    directories["collocation_indices"] = \
        os.path.join(directories["collocation_files"],
                     "2007_Index_CALIPSO-MODIS-CloudSat")

    directories["CALIPSO_01km_data"] = \
        os.path.join(directories["group_data"],
                     "2007_Data_CALIPSO_01kmCLay")
        
#    directories["CALIPSO_05km_data"] = \
#        os.path.join(directories["group_data"],
#                     "2007_Data_CALIPSO_05kmCLay")
        
    del directories["group_data"]
    del directories["collocation_files"]
        
    return directories
