#!/bin/bash

# This function writes a script for driving HEG to merge LAI HDF files 
# and convert them to GeoTIFF files.

DIR_NAME=$(cd "$(dirname "$0")";pwd)

echo $DIR_NAME

files1104=$(ls ${DIR_NAME}/heg.h11v04/)
files1105=$(ls ${DIR_NAME}/heg.h11v05/)
files1204=$(ls ${DIR_NAME}/heg.h12v04/)
files1205=$(ls ${DIR_NAME}/heg.h12v05/)

ARR1=($files1104)
ARR2=($files1105)
ARR3=($files1204)
ARR4=($files1205)

declare -a ARR1
declare -a ARR2
declare -a ARR3
declare -a ARR4

for ((i=0; i<=$[${#ARR1[*]}-1]; i ++))
do
    filename1=${ARR1[${i}]}
    filename2=${ARR2[${i}]}
    filename3=${ARR3[${i}]}
    filename4=${ARR4[${i}]}
    length=${#filename1}
    outname=${filename1:0:17}

    cat << _EOF_ > MyParameter.prm

NUM_RUNS = 1

BEGIN
NUMBER_INPUTFILES = 4
INPUT_FILENAMES = ${DIR_NAME}/heg.h11v04/${filename1}|${DIR_NAME}/heg.h11v05/${filename2}|${DIR_NAME}/heg.h12v04/${filename3}|${DIR_NAME}/heg.h12v05/${filename4}
OBJECT_NAME = MOD_Grid_MOD15A2H|
FIELD_NAME = Lai_500m|
BAND_NUMBER = 1
SPATIAL_SUBSET_UL_CORNER = ( 41.0 -80.0 )
SPATIAL_SUBSET_LR_CORNER = ( 37.0 -76.0 )
OUTPUT_OBJECT_NAME = MOD_Grid_MOD15A2H|
OUTGRID_X_PIXELSIZE = 0.0041666733029590754
OUTGRID_Y_PIXELSIZE = 0.0041666733029590754
RESAMPLING_TYPE = NN
OUTPUT_PROJECTION_TYPE = GEO
ELLIPSOID_CODE = WGS84
OUTPUT_PROJECTION_PARAMETERS = ( 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  )
OUTPUT_FILENAME = ${DIR_NAME}/tif_files/${outname}.tif
SAVE_STITCHED_FILE = NO
OUTPUT_STITCHED_FILENAME = ${DIR_NAME}/tif_files/${outname}.hdf
OUTPUT_TYPE = GEO
END
_EOF_
    subset_stitch_grid -p MyParameter.prm
    rm -f MyParameter.prm
done

