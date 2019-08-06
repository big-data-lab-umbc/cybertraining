Team 4 Project of the CyberTraining program at UMBC in 2019 (http://cybertraining.umbc.edu/)

**Title**: Deep Learning Based Mineral Dust Detection and Feature Selection

**Team members**: Ping Hou, Peng Wu

**Mentors**: Aryya Gangopadhyay

**Abstract**: Dust storm affects human health and the environment. In this study, we develop deep learning models to identify dust from cloud and surface using MODIS observations and CALIPSO data. We also identified the best subset of channels for dust detection by a shuffling procedure and a genetic algorithm. Results show the important features determined by the two methods are very similar. And the genetic algorithm selected a subset of features that have a comparable performance with the model with all features, proving the effectiveness of the genetic algorithm. The chosen subset will reduce future data collection efforts for dust detection.

**Structures of implementation**:

"Deeplearning_dust_detection_march2014.ipynb" contains the deep learning model we developed based on CALIPSO data in march 2014. It also has the code for the shuffling procedure to identify the importance of features. 

"Deeplearning_dust_detection_feature_selection_parallel.py" and "run_parallel.slurm" are the codes for feature selection using a genetic algorithm to select a subset of features for duct detection. 

**Instructions on how to run the code**

Users can run "Deeplearning_dust_detection_march2014.ipynb" on Jupyter notebook. 

Users can run "Deeplearning_dust_detection_feature_selection_parallel.py" by "run_parallel.slurm" on Taki. Users can change the parameters (e.g.,the population size, the number of generations) at the begnining of the code.  
