* 1D/3D Radiative Transfer Simulation using SHDOM
Meng Gao <mgao@umbc.edu>
08/03/2018

The work is from an collaboration effort with
Scott Hottovy <hottovy@usna.edu>
Cui, Yunwei" <ycui@towson.edu>
Zhibo Zhang <zzbatmos@umbc.edu>
Chamara Rajapakshe <charaj1@umbc.edu>

* NOTE:
The following discriptions involve multiple steps. Each step includes
the definition of a template input file, python and bash scripts to customize
the input file, and scripts to submit the jobs. Eventually after all the data are 
computed, the finale resulst will be extracted from the outputfile, 
and analyzed using python cripts. All the necessary functions, 
and example data processing procedures are summarized in the ipython notebook. 

Please update the data path and make sure the scripts can locate your data. 

* SHDOM code installation 
1. SHDOM source code can be downloaded from:
http://nit.colorado.edu/shdom.html
2. Please follow the instruction provided with the source code for installation.
3. The current study is based on the unpolarized version of SHDOM.

* Single scattering property generation
Scripts and sample results to generation single scattering properties are in /shdom_run/scat_table
please make the results locatable for the radiative transfer simulation. 

* 1D and 3D RT simulation
1. All the sripts, template input file, script for job submission are provided in /shdom_run/cloud_2d_modis4
2. Please refer to SHDOM manual for the detailed information on the SHDOM input parameters. 
3. After finish the run, please copy the output file into a data folder for the following data analyzing procedure. 

* Data analyzing
1. A ipython notebook is provided in /shdom_analyze.
2. The notebook demonstrates my procedure to analyze the data with relevant functions provided. 
3. Following the steps in the notebook, the usage of the functions should be self-explanatory.
4. Functions to plot the 1D and 3D reflectance, the CER/COT look up table, and retrieval results are also provided. 
5. One set of the test results with the plots are provided as examples in the notebook.
