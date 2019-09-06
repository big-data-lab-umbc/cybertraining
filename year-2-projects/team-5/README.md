Team 5 Project of the CyberTraining program at UMBC in 2019 (http://cybertraining.umbc.edu/)

**Title**:
DUST DETECTION IN SATELLITE DATA USING CONVOLUTIONAL NEURAL NETWORKS (CNN)

**Team members**: 
Changjie Cai, Jangho Lee, Yingxi Rona Shi, Camille Zerfas

**TA**: 
Pei Guo

**Mentors**: 
Zhibo Zhang

**Abstract**:
Atmospheric dust is known to cause health ailments and impacts earth’s climate
and weather patterns. Due to the many issues atmospheric dust contributes to, it is
important to study dust patterns and how it enters the atmosphere. In the past, many
scientists have used satellite data and physical-based algorithms to detect and track
dust, but these algorithms have many shortcomings. Herein, we consider Convolutional
Neural Networks to classify dust in satellite images to try to improve the accuracy of
dust detection. We describe the satellite data used, discuss the model structures, and
provide results for the models built. These models show promising preliminary results.

**Structures of implementation**:
python_01_datapreprocessing:  
  01a_read_and_merge_data.py ~ 01d_read_and_merge_data.py: script files reads the raw VIIRS and CALIPSO files and extract training / testing dataset with 5x5 moving window. Each of a~d scipt files reads raw data by season.   
  02_normalize.py reads extracted data from previous script code and normalize them.  
  03_count_classes.py counts the number of classes in the data.

python_02_model:
  Each of the script files in this directory is a CNN model.
 
 python_03_result:
  Scripts in this directory is used for analyzing the result of each models
  
**Instructions on how to run the code**
The codes should be ran swquentially with numbers. Data is not contained in this repository, due to its large set of data.
