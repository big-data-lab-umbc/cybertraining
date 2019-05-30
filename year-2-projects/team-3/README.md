Team 3 Project of the CyberTraining program at UMBC in 2019 http://cybertraining.umbc.edu/
#### Title: An Approach to Tuning Hyperparameters in Parallel - A Performance Study 
#### Team members: Charlie Becker; Bin Wang; Will Mayfield; Sarah Murphy
#### Mentors: Dr. Matthias Gobbert; Carlos Barajas


This is a working example of a performance study completed on the 'Taki' HPC cluster at UMBC. It uses a combination of popular Python modules for hyperparameter tuning in parallel. The data and base model configuration is borrowed from the Machine Learning in Python forEnvironmental Science Problems AMS Short Course, provided by David John Gagne from the National Center for Atmospheric Research.  The repository for that course can be found at [https://github.com/djgagne/ams-ml-python-course]

Below, are some brief directions to reproduce the results in the technical report. Full results are seen and discussed in [Technical\_Report.pdf](https://github.com/big-data-lab-umbc/cybertraining/blob/master/year-2-projects/team-3/Technical_Report.pdf) 
### Workflow
After cloning this directory, first run data\_download.py which will download the data into a data direcory.

Next, run preprocess.py which will preproces the data and augement it to give a balanced dataset.  It will create and place .npy files into the data directory for easy access.

Next, you can run either submit\_2013.slurm (2013 partition), submit\_2018.slurm (2018 partition) or submit\_gpu.slurm (2018 GPU nodes) to submit the performance study across the cluster.  These scripts call run\_2013.py, run\_2018.py and run\_gpu.py respectively, which is where additional SLURM argumentes are defined, such as the number of nodes and hyperparameters.  Specifically, cluster.scale(x) will refer to the number of nodes desired.

Output for the study will be delivered to slurm-2013.out, slurm-2018.out or slurm-gpu.out with error logs being delivered to slurm-xxxx.err. Additionally, each process within each node will prodeuce training output in slurm-jobID.out, though this probably won't be useful.
### Data augmentation
RandomOverSampler class from imblearn.over\_sampling was used to oversample the minority classes (non-tornadic data) fed into the deep neural network. This is done to achieve an approximate 50/50 class split within the training data; which began as approximately a 95/5 split. The relevant script is dnn.py

For convolutional neural network, the input data are tensor images. We augment the minority classes (non-tornadic images) by duplicating, shuffling, and transforming the images through small angle rotation but keeping the labels unchanged. This can be done in real time while training the model via ImageDataGenerator from Keras or at the preprocessing stage using skimage.transform.rotate before data are feeding into the model. In the code, this can be selected via two parameters 'augmentation' and 'on\_the\_fly'. For example, if augmentation==True, and on\_the\_fly==False, this means that the augmented data is generated before training. The working script is cnn.py

Overall, we did not see a significant spike in performance when using transformed data as opposed to resampled data only.  However, augmentation is highly specific to the dataset and will have varying benefits dependent on each specfic dataset.   
