Team 4 Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

# Use of Deep Learning to Classify Compton Camera BasedPrompt Gamma Imaging for Proton Radiotherapy

# Team Members
* Jonathan N. Basalyga 
* Gerson C. Kroiz

# Mentors
* Carlos A. Barajas, M.S.
* Matthias Gobbert, Ph.D.

# Clients
* Paul Maggi, Ph.D.
* Jerimy Polf, Ph.D.

# Abstract
Real-time imaging has potential to greatly increase the effectiveness of proton
beam therapy for  cancer  treatment.   One  promising  method  of  real-time
imaging  is  the  use  of  a  Compton camera to detect prompt gamma rays, which
are emitted along the path of the beam, in order to reconstruct their origin.
However,  because of limitations in the Compton camera's ability to  detect
prompt  gammas,  the  data  are  often  ambiguous,  making  reconstructions
based  on them unusable for practical purposes.  Deep learning's ability to
detect subtleties in data that traditional models do not use make it one possible
candidate for the improvement of classification of Compton camera data.  We show
that a suitably designed neural network can reduce false detections and
misorderings of interactions, thereby improving reconstruction quality.

# Instructions on how to run the code
For instructions on how to use the code from start to finish please see
[the python notebook](./howToPGML.ipynb)

Note that there is no guidance on how to train the network.
All hyperparameters for training are in [params.json](./params.json).
1. Download and load the h5 file with Keras
2. Reset the weights using the advice stated on [their github](https://github.com/keras-team/keras/issues/341)
3. Use the hyperparameters from [params.json](./params.json)
4. Then call `model.fit(inputData, outputData)`


