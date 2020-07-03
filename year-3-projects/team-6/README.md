Team 6 Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

**Title**: Studying Arctic Sea Ice Retreat Using Data-Driven Causality Discovery Approaches

**Team members**: Yiyi Huang, Matthäus Kleindessner, Alexey Munishkin, Debvrat Varshney

**Mentors**: Jianwu Wang, Pei Guo

**Abstract**: 

**Structures of implementation**:

**Instructions on how to run the code**:

For plotting TCDF, NOTEARS and DAG-GNN graphs, you need the [igraph](https://igraph.org/python/) package.

###### TCDF

In order to run the code, you need to download the TCDF code from [this link](https://github.com/M-Nauta/TCDF). The instructions to run the code is given on that link, but below are specific instructions for this project.

We run the code on [taki high performance computing faculty](https://hpcf.umbc.edu/) so first you need to
```
module load networkx
module --ignore-cache load "Python/3.7.6-intel-2019a"
```
then you can run the TCDF code
```
python runTCDF.py --data report_data/combined_decomposed_drop_temp_all_norm_1980_2018.csv --kernel_size k# --hidden_layers h#
```
where `report_data/combined_...` is the same csv data used for the NOTEARS and DAG_GNN code. Note: `k#` and `h#` are used for the hyperparameter sensitivity study where for this project we used `k# = {2,4,6}` and `h# = {0,1,2}`. Note: replace `k#` and `h#` with numeric integer values.

After you run the TCDF code you will get a text display, e.g.
```
python runTCDF.py --data report_data/combined_decomposed_drop_temp_all_norm_1980_2018.csv --kernel_size 2 --hidden_layers 0
...
...
===================Results for combined_decomposed_drop_temp_all_norm_1980_2018.csv ==================================
v10m causes u10m with a delay of 0 time steps.
u10m causes v10m with a delay of 0 time steps.
========================================================================
```
which can (manually) converted to a numpy adjacency matrix
```
([[0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,1,0,0,0,0],
 [0,0,0,0,0,0,1,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,0,0,0,0,0,0,0,0,0]])
```
where the going from left to right or up to down is `['HFLX','SW','LW','SLP','Precip','RH','u10m','v10m','sea_ice','CC','CW','GH']`, which are aliases for atmospheric and sea ice variables (this is also used by NOTEARS and DAG-GNN)

* longwave: LW
* shortwave: SW
* tot_precip: Precip
* u10m: u10m
* v10m: v10m
* sea_ice: sea_ice
* GH_mean: GH
* SLP: SLP
* RH: RH
* cloud_cover: CC
* cloud_water: CW
* heat_flux: HFLX

now we can plot the data (though you have to manually adjust the numpy matrices in the code)
```
python plot_graph_for_TCDF.py
```
and for the hyperparameter sensitivity study (make sure to manually adjust the matrices here too)
```
python compute_normHamming_for_TCDF.py
```

###### NOTEARS

In order to run the code, you need to download the NOTEARS code from https://github.com/xunzheng/notears. We used Version 2.1 of the NOTEARS code.

Then just run 
```
run_notears_on_data_STATIC.py 
```
or 
```
run_notears_on_data_TEMPORAL.py 
```

###### DAG-GNN

Copy the `src` folder from [this link](https://github.com/big-data-lab-umbc/DAG-GNN) and place it inside the DAG-GNN folder here. Then put the detrended and deaseasonalized datasets inside the `data` folder in csv format. Your DAG-GNN folder structure should look like this -
* DAG-GNN
  * data
    * file1.csv
    * file2.csv
    * .
    * .
  * plots
    * plot_graph.py
    * plotGNN.py
  * src
    * train.py
    * utils.py
    * modules.py

To run DAG-GNN, just run the following command from inside the `src` folder
```
python train.py --filename=file1.csv --epochs=50 
```
Your DAG will be generated at this location: `src/<filename>__epochs<no. of epochs>/predG`. For example, for the command above, your DAG will be a file called "predG" generated in `src/file1__epochs50/`.

Further, DAG-GNN can be customized with a whole list of arguments, but the `filename` argument is compulsory. If the number of epochs is not defined, then its default value will be considered, which is 200.

To plot the graphs, simply run the following command from inside the `plots` folder
```
python plotGNN.py 
```
Each graph will be generated in the folder containing the corresponding "predG".
