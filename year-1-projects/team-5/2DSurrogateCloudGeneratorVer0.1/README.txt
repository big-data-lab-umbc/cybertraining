# Project Title

Surrogate2DCloudGenerator.m

Two dimensional surrogate cloud model. 
This code generates a 2 dimensional (lattitude and longitude) surrogate cloud based off given MODIS data of
 cloud optical thickness. Then the variables of cloud effective radius (CER) and cloud top height (CTH) are 
derived from COT using atmospheric physics relations. See technical report. 

## Getting Started

The code is written in matlab as a script. To run

### Prerequisites

MATLAB
-Surrogate2DCloudGenerator.m
-periodogram2.m
-MODIS data files of COT and CER. Example is given in the directory MODISCOT. 

### Installing

No installation necessary. 

## Running the tests

Open the Surrogate2DCloudGenerator.m file and hit the "Run" button. 

## Authors

Scott Hottovy (https://www.usna.edu/Users/math/hottovy/index.php)

## Acknowledgments

* Collaborators Yunwei Cui and Meng Gao
* Code is built from a description form the paper Hogan & Illingworth (1999)
* 2D power spectrum code used from Wang Xianju. 
