# K-means_Clustering4CloudHistogram
K-means clustering method applied to Cloud 2D joint histogram

---
### 1. Preparing input data
MODIS_raw-hdf2binary.py3.py
: Transforming raw hdf file format MODIS data to grided cloud 2D joint histogram binary format
: Output: [nt(days), 180(lats), 360(lons), 7(CTP), 6(COT)], 4-byte float (C- or Python-structure, so COT first, and time last)

MODIS_data_extract.py3.py
: Sub-sampling of binary file by selecting time range and/or regional range

### 2.1 Running_K-means_Ver_Fortran+CSH
kmeans_omp_compiling_v3a_dpout_TR_b42.csh
: This is ready to run script if input file is ready and set the file location.
: Before run the code, prepare a directory for output, and check the location of output directory in the c-shell script
: This c-shell script embed fortran program using OpenMP. 
: The fortran program is tested with Intel fortran and gfortran compliers

### 2.2 Running_K-means_Ver_Python+F2py
Under the work

### 3 Check_output
cent_init_display_42bins_obs.py3.py
: This is a python3 program to display centroid (upto k=8). 

Additional programs like checking geographical distribution of "Relative Frequency of Occurrence (RFO)" are under work. 

