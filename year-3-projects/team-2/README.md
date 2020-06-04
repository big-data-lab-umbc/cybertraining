Team 2
======

Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/).


Title:
------

***Studying Anomalous Discrepancies between MODIS and CALIOP Cloud Observations***


Team members:
-------------

 * Christine Abraham<sup>2</sup>
 * Olivia Norman<sup>1,3</sup>
 * [Erick Shepherd](https://github.com/ErickShepherd)<sup>1,3</sup>


Mentors:
--------

 * **Faculty advisor:** [Dr. Zhibo Zhang](https://www.researchgate.net/profile/Zhibo_Zhang2)<sup>1,3</sup>
 * **Research assistant:** [Jianyu Zheng](https://www.researchgate.net/profile/Jianyu_Zheng3)<sup>1,3</sup>


Affiliations:
-------------

**[University of Maryland, Baltimore County (UMBC)](https://www.umbc.edu/):**
 1. [Aerosol, Cloud, Radiation-Observation and Simulation (ACROS) Group](https://acros.umbc.edu/), 
 2. [Department of Information Systems](https://informationsystems.umbc.edu/)
 3. [Department of Physics](https://physics.umbc.edu/)


Abstract:
---------

When examining collocated data from the A-Train satellite constellation, there
are a notable number of clouds that CALIOP identifies as transparent but which
MODIS paradoxically reports have a high cloud optical thickness (COT). Our team
is investigating two hypotheses in an effort to explain the occurrence of these
anomalies: that they could be MODIS COT retrieval errors due to the
misclassification of high albedo surfaces, such as snow and sea ice, as clouds;
or that they could be clouds which are misidentified as having a high COT due
to 3D radiative effects. In investigating the former hypothesis, we collocated
the single-layer over water anomalies with NSIDC AMSR-E sea ice observations
using a k-nearest neighbors (k-NN) algorithm and determined that around 50% of
such anomalies occur over areas with high sea ice concentrations (95-100%).
With half of the results showing no correlation with high albedo surfaces,
further research is required to account for the difference. We have taken
preliminary steps toward exploring whether the cloud 3D radiative effects
hypothesis might explain the remaining anomalies.


Structures of implementation:
-----------------------------

 1. **Identify the MODIS-CALIOP collocation anomalies:**
 
    1. Load the collocation data.
    
    2. Mask the single-layer clouds.
    
    3. Mask the subset of clouds by underlying CERES surface (water or land).
    
    4. Mask the subset of clouds by phase (water vapor or ice).
    
    5. Mask the subset of clouds which CALIOP flags as transparent.
    
    6. Mask the subset of clouds which MODIS reports have a COT > 150.
    
    
 2. **Identify correlations:**
 
    1. Generate a scatter plot map of where the anomalies occur.
    
    2. Set the color of the scatter points to the solar zenith angle.
    
    3. Generate an intensity plot of transparent cloud fraction vs. COD and SZA.
    
    
 3. **Explore the sea ice hypothesis:**
 
    1. Collocate the anomalies with NSIDC sea ice concentration data.
    
    2. Plot a histogram of anomaly counts vs. sea ice concentration.
    
    
 4. **Explore the 3D radiative effects hypothesis:**
 
    1. Extract the MODIS cloud top height data.
    
    2. Derive the CALIOP cloud top height.
    
    3. Compute the slopes of the cloud top heights about the anomalies.
    
    4. Determine whether the anomaly occurs on the illuminated or shaded side of
       the cloud.


Instructions on how to run the code:
------------------------------------

At the stage of development seen in this snapshot of our work, our codebase
consists of a series of separate scripts and utility modules, with no singular
interface or main entry point. As is, scripts are intended to be run separately
for the singular purposes for which they were designed. Our workflow was to
design working prototype code to perform the needed tasks first and refactor it
later.

**For the most recent version and documentation of our software, visit [our primary GitHub repository](https://github.com/ErickShepherd/modis_caliop_anomaly_analysis):**
 * https://github.com/ErickShepherd/modis_caliop_anomaly_analysis
