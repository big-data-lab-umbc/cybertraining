Team 2
======
Project of the CyberTraining program at UMBC in 2020 (http://cybertraining.umbc.edu/)

Title:
------

**Machine Learning for Retrieving Cloud Optical Thickness from Observed Reflectance: 3D Effects**

Team members:
------------
 
 * Kallista Angeloff
 * Kirana Bergstrom
 * Tianhao Le
 * Chengtao Xu

Mentors:
------- 
 * Prof. Zhibo Zhang

Abstract: 
--------
Clouds are inherently 3D, and simulating radiative transfer (RT) properties accurately requires models that take their 3D effects into account. Because 3D models are complex and computationally expensive, RT models often use simplified 1D models to retrieve cloud properties, which suffer from inaccuracies due to 3D effects. Recent advancements in deep learning may lead to a retrieval algorithm that is capable of taking these effects into account. We will develop a deep-learning based cloud property retrieval algorithm that is able to reconstruct the 3D structure of clouds based on observed cloud radiative signatures.

Structures of implementation:
----------------------------
1. shdom 
   * 3D RT code

2. Frac_cloud
   * frac_cloud.py

      source code for generate fractional cloud
   * profile.py

      generate profiles for RT
   * RT_wrap_parallel.py

      Run RT 
   * combine.py

      Combine all RT outputs for ML input

3. DATA

   direcotory for ML inputs

4. ML 
   * cnn.py
     serial version
   * par_cnn.py
     parallel version
 
