# (ECO) Evolutionary clustering of moving objects
## Introduction
The folder "code" holds the source code of our method.
The folder "data" holds a sample of CD dataset used in our paper.

## Environment Preparation
  - Vistual Studio 2019
  - Windows 10 Pro, 64-bit operating system
  - Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz 
 
## Datasets Description
  - File CD_example.txt gives the input format of ECO.
  - Each item of tha data (each row) has 16 properties (1-16).
  - We use properties 2, 3, 4, and 15, which are the ID of a moving object, longitude of its location, lantitude of its location and its arriving timestamp .
  - Data from each of two time steps is separeted by "= = = = = = = = = = = = = = = =".
  
## Code Description
  - Header.h stores the settings of parameters and the definitions of data structrues used in the code.
  - get_input.h claims all the functions of the code.
  - get_input.cpp defines all the functions of the code.

## Usage 
  - Unzip "ECO.zip"
  - Create a project using files in the folder "ECO" and run 
 Note: the route of the dataset is in evolutionary_clustering.cpp

## Parameter Setting
  - In Header.h, 
    - "delta" and "rho" are the parameters "\delta" and "\rho" for generating minimal groups
    - "xi" is the parameter "minPts" in DBSCAN
    - "alpha" is the parameter "\alpha" in evolutionary clustering
    - "mu" is the speed constraint of moving objects
    - "delta_epsilon"  is the parameter "\Delta_\varepsilon" for adjusting "\varepsilon" at each time step
  - In evolutionary_clustering.cpp, 
    - "epsilon"  is the initial value of the parameter "\varepsilon" in DBSCAN
  -Default setting
    - delta=1000m, rho=5, xi=3, alpha=0.1, mu=11.1m/s, delta_epsilon=50m and epsilon=1200m
