Team 7 Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

**Title**: Multi-sensor dust detection using unsupervised machine learning

**Team members**: 
Julie Bessac 	(Mathematics and Computer Science Division, Argonne National Laboratory), 
Ling Xu 		(Department of Mathematics, North Carolina A&T State University), 
Manzhu Yu 	(Department of Geography, Penn State University)

**Mentors**: 
Dr. Aryya Gangopadhyay (Department of Information Systems, University of Maryland  Baltimore County)

**Abstract**: 
Dust and sandstorms originating from Earth’s major arid and semi-arid desert areas can significantly affect the climate system and health. Many existing methods use heuristic rules to classify on a pixel-level regarding dust or dust-free, but these heuristic rules are limited in applicability when the study area or the study period is changed. Based on the multi-sensor collocation dataset generated from CyberTraining 2019, we sought to utilize unsupervised machine learning techniques to detect and segment dust in multispectral satellite imagery. 

**Structures of implementation**:

1.download_viirs_samples.py: use collocated calipso-viirs dataset to download VIIRS granule based on filename

Experiments using 256 by 256 pixels subsets of VIIRS imagery collocating with CALIPSO tracks:

2.cal_pixel_label_256.py: generate 256 by 256 images/data: mask, predictor, figure, composite, and landtype

3.kmeans_one_image.py: use K-means clustering to cluster one image within the test set

4.kmeans_average_accuracy.py: use kmeans clustering method to cluster all images within the test set and calculate average accuracy

Experiments using the entire subsets of VIIRS imagery collocating with CALIPSO tracks:

5.test_entire_granule.py: experiment test on a VIIRS-CALIPSO granule subset from extracting data, K-means, accuracy evaluation in one python code

**Instructions on how to run the code**

All codes able to run on taki
Code location: /umbc/xfs1/cybertrn/cybertraining2020/team7/research/VIIRS-SIPS/
