Team 8 Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

**Title**: Tornado Prediction using Environmental Sounding Data: Comparing Random Forest to CNN

**Team members**: Brice Coffer (Department of Marine, Earth, and Atmospheric Science, North Carolina State University), Michaela Kubacki (Department of Mathematics, Middlebury College), Yixin Wen (Cooperative Institute for Mesoscale Meteorological Studies, University of Oklahoma, and NOAA/National Severe Storms Laboratory, Norman, Oklahoma), Ting Zhang (Department of Mathematics and Computer Science, McDaniel College)

**Mentors**: Matthias Gobbert (Department of Mathematics and Statistics, University of Maryland, Baltimore County) and Carlos Barajas (Department of Mathematics and Statistics, University of Maryland, Baltimore County)

**Abstract**: Tornadoes pose a forecast challenge to National Weather Service forecasters because of their quick development and potential for life-threatening damage. The use of machine learning in severe weather forecasting has recently garnered interest, with current efforts mainly utilizing ground weather radar observations. In this study, we investigate machine learning techniques to discriminate between nontornadic and tornadic storms solely relying on the Rapid Update Cycle (RUC) sounding data that represent the pre-storm atmospheric conditions. This approach aims to provide for early warnings of tornadic storms, before they form and are detectable by weather radar observations. Two machine learning methods tested in our project are Random Forest (RF) and Convolutional Neural Network (CNN). Performance testing of RF using various ranges of hyperparameters results in an overall accuracy score of 70.14%, but the accuracy of significantly tornadic class prediction is only 23.84%. The CNN model results in an overall accuracy score of 67.84%, but the accuracy for significantly tornadic storms is only 26.69%. The higher accuracy in the RF and CNN models for the majority class of nontornadic supercells suggests that the imbalanced dataset is a meaningful contributor to the lower accuracy for tornadic storms. After applying the simple method of randomly undersampling (oversampling) the majority (minority) class, the accuracies of significantly tornadic class prediction of RF and CNN are enhanced to 65.85% and 36.01%, respectively.

**Structures of implementation**:
1. CNN/cnn1D.py: creates a CNN model and evaluate it
2. CNN/convert_variable.py: converts base state variables in RUC soundings to other possible representations
3. CNN/create_model.py: creates two complicated CNN models
4. CNN/evaluate_model.py: evaluates CNN models and plots charts to show accuaracies over epochs
5. CNN/imbalanced_data.py: handles imbalanced datasets issue using undersampling and oversampling approaches
6. CNN/interpolateRUCsounding.py: interploates feature values with varying height levels over a fixed vertical height grid above ground level (AGL)
7. CNN/normalize_variables.py: normalizes each variable from 0 to 1
8. CNN/readRUCsounding.py: uses to read in RUC sounding and returns a dictionary of metadata and an 3D array of sounding data, structured as a 20194x37x7 numpy array
9. RandomForest/RF_model.py: creates, trains, and tests a Random Forest classifier for the RUC sounding data for supercell storms that has been interpolated over 37 fixed heights
10. RandomForest/interpolate_avg_height.py: converts data with varying height levels into data with a fixed height grid comprised of the average heights at each of the 37 levels and interploates feature values over the fixed height grid
11. RandomForest/RF_oversample.py: creates and tests a Random Forest model that handles imbalance datasets issue using oversampling approach
12. RandomForest/RF_undersample.py: creates and tests a Random Forest model that handles imbalance datasets issue using undersampling approach

**Instructions on how to run the code**
All codes able to run on taki
Code location: /umbc/xfs1/cybertrn/cybertraining2020/team8/research
