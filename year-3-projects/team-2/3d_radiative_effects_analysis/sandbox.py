from pyhdf.SD import SD, SDC
import vertical_feature_mask as vfm

filename = ("C://Users//Erick Shepherd//Desktop//NSF REU//Data//"
            "2010_Data_CALIPSO_01kmCLay//CAL_LID_L2_01kmCLay-ValStage1-V3-01."
            "2010-01-01T00-22-28ZD.hdf")
file     = SD(filename, SDC.READ)
data     = file.select("Feature_Classification_Flags").get()
idf_data = vfm.extract_features(data)
