# GenBinDat_wCorr
Generates random binary data with a chosen values of correlation between the variables. 

The correlation of the generated datasets is centered around the specified correlation with a standard deviation of the order of 1/sqrt(N), where N is the number of datapoints of the generated datasets.

## Requirements

The code uses the C++11 version of C++.

It also uses the two header-only C++ libraries: Eigen and LBFGS++ (L-BFGS solver).
Remember to specify the paths to the two libraries to your compiler (see `-I` option below).

**To compile:**  `g++ -I /usr/local/include/eigen3/ -I ./LBFGSpp-master/include -O2 main.cpp Models.cpp ReadDataFile.cpp BoltzmannLearning.cpp FitFromCorrMatrix.cpp Generate_data_exact.cpp`

**To execute:** `./a.out`


