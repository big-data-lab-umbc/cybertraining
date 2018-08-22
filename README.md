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
: This is very similar code to mimic previous Fortran+CSH version.  
: There are two modules: 1. python class (k_means_class_py3.py), 2. Fortran Module (k-means_mod.f90)  
: The fortran module should be compiled first, and shell script file contains f2py command and options (f2py_k-means_module_omp.csh)  
: An example code to run K-means clustering is given (run_k-means_py3.py) 

### 2.3 MPI\_Cython
This code uses MPI and OpenMP to squeeze as much performance as possible from a system. 
It is assumed that you already have a functioning MPI distribution and a compatible C compiler.
Simply adjust and run `./build_k_means.sh` to create the shared library for python.
To run:

1. Adjust `indir` and `infile` accordingly in `run_k-means_py3.py`
2. Adjust the `kmeans.slurm` sample file to match your clusters needs
3. Submit the file to your scheduler.
4. If you aren't using an NFS mount or something similar you'll need to ensure that the following 
files are transferred to the cluster: `shared library`, `run_k-means_py3.py`, `k_means_class_py3.py`.
    5. The file output by this can be post-processed by `3`.

### 2.4 Run K-Means Using SparkML Library 
This is new implementation using SparkMLlibrary.

How to run it:
    1. First run test2read_clustering_input_py3.py and generate a .csv file as Spark programs' input. 
:
    2. Submit spark_ml_kmeans.py to spark. spark-submit spark_ml_kmeans.py seedid. (seed id is the kmeans seed id for repeatable reason)
    
:  e.g. spark-submit spark_ml_kmeans.py 1

:   3. The output is printed in the console.
 
### 3 Check_output
cent_init_display_42bins_obs.py3.py  
: This is a python3 program to display centroid (upto k=8).   
cent_init_display_42bins_obs.k9+.py3.py  
: This is a python3 program to display centroid (9<k<=12). 

~~Additional programs like checking geographical distribution of "Relative Frequency of Occurrence (RFO)" are under work.~~

CR_num_map_class_py3.py   
run_CR_num_map_class.py3.py   
: Class and running script to get CR_num file.  
: For every data points in input data, determine the closest centroid (called as assigning), and return the cluster number  
: Hence, if input data format=[nday,nlat,nlon,nelem], output format=[nday,nlat,nlon], numpy.int16  
: -1: missing, 0: clear sky (no clouds), 1 to k: cluster number.  

rfo_map_display_class_py3.py  
run_rfo_map_display.py3.py  
: Class and running script to draw RFO map from CR_num file  

check_rfo_map_similarity.py3.py  
: Shows similarity between cluster rfo maps and betwen centroids, in terms of correlation and RMSD, using "rfo_map_display_class"
 

