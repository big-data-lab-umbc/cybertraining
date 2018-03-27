# K-means_Clustering4CloudHistogram
K-means clustering method applied to Cloud 2D joint histogram

---
## 1. Version_Fortran+Csh
This is the original version using C-shell and Fortran program. 
  
### 1.1 Preparing input data
Currently, the python3 program extract part of cloud histogram variable for 1 year from "pre-processed binary" file. The program to make an input data from original hdf source is under work.

### 1.2 Running_K-means
This is ready to run if input file is ready and set the file location.
Before run the code, prepare a directory for output, and check the location of output directory in the c-shell script

### 1.3 Check_output
Currently, a python3 program to display centroid (upto k=8) is given. Additional programs like checking geographical distribution of "Relative Frequency of Occurrence (RFO)" are under work. 

---
## 2. Version_Python_with_F2py
Currently working on. 
