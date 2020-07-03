Team 6 Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

**Title**: Studying Arctic Sea Ice Retreat Using Data-Driven Causality Discovery Approaches

**Team members**: Yiyi Huang, Matthäus Kleindessner, Alexey Munishkin, Debvrat Varshney

**Mentors**: Jianwu Wang, Pei Guo

**Abstract**: 

**Structures of implementation**:

**Instructions on how to run the code**:

For plotting NOTEARS and DAG-GNN graphs, you need the [igraph](https://igraph.org/python/) package.
###### TCDF

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

Copy the `src` folder from [this link](https://github.com/big-data-lab-umbc/DAG-GNN) and place it inside the DAG-GNN folder here, in our team's repo. Then keep the detrended & deaseasonalized dataset inside a `data` folder in csv format. Your DAG-GNN folder structure should look like this -
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

To run DAG-GNN, just run
```
python src/train.py --filename=file1.csv --epochs=50 
```
Your DAG will be generated at this location: `src/<filename>__epochs<no. of epochs>/predG`. For example, for the command above, your DAG will be a file called "predG" generated in `src/file1__epochs50/`.

Further, DAG-GNN can be customized with a whole list of arguments, but the `filename` argument is compulsory. If the number of epochs is not defined, then its default value will be considered, which is 200.

To plot the graphs, simply run
```
python plots/plotGNN.py 
```
Each graph will be generated in the folder containing the corresponding "predG".
