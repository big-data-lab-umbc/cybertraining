#!/bin/bash

# This function writes a HEG Parameter script 
# for converting HDF regional files (h12v04) to GeoTIFF files

DIR_NAME=$(cd "$(dirname "$0")";pwd)

echo $DIR_NAME

files=$(ls ${DIR_NAME}/heg.h12v04/)

for filename in $files
do
    length=${#filename}
    outname=${filename:0:${length}-4}
    cat << _EOF_ > MyParameter.prm

NUM_RUNS = 1

BEGIN
INPUT_FILENAME = ${DIR_NAME}/heg.h12v04/${filename}
OBJECT_NAME = MOD_Grid_BRDF|
FIELD_NAME = Albedo_BSA_shortwave
BAND_NUMBER = 1
SPATIAL_SUBSET_UL_CORNER = ( 41.000000000 -79.5 )
SPATIAL_SUBSET_LR_CORNER = ( 39.999999996 -76.0 )
RESAMPLING_TYPE = BI
OUTPUT_PROJECTION_TYPE = GEO
ELLIPSOID_CODE = WGS84
OUTPUT_PROJECTION_PARAMETERS = ( 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  )
OUTPUT_PIXEL_SIZE = 15.000023890652672
OUTPUT_FILENAME = ${DIR_NAME}/tif_files/heg.h12v04/${outname}.tif
OUTPUT_TYPE = GEO
END

_EOF_
    resample -p MyParameter.prm
    rm -f MyParameter.prm
done

