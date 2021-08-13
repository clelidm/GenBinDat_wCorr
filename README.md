# GenBinDat_wCorr
This program generates random binary data with a chosen values of correlation between the variables. To do so it fits a fully connected pairwise spin models to reproduce the specified correlation patterns, and then generates data from the fitted model. In the dataset, datapoints are independently sampled from the same multivariate probability distribution (i.e., from the model).

The program can also generates data from any specified pairwise models fitted on a given dataset.

The correlation of the generated datasets is centered around the specified correlations (or the correlations of the original dataset) with a standard deviation of the order of 1/sqrt(N), where N is the number of datapoints of the generated datasets.

## Requirements

The code uses the C++11 version of C++.

It also uses the two header-only C++ libraries: Eigen and LBFGS++ (L-BFGS solver).
Remember to specify the paths to the two libraries to your compiler (see option `-I` below).

**To compile:**  `g++ -I /usr/local/include/eigen3/ -I ./LBFGSpp-master/include -O2 main.cpp Models.cpp ReadDataFile.cpp BoltzmannLearning.cpp FitFromCorrMatrix.cpp Generate_data_exact.cpp`

**To execute:** `./a.out`

## Examples

All the functions that can be used from `int main()` are declared at the beginning of the `main.cpp` file and described below.

For hands-on and simple tests of the program, check the examples in the `int main()` function of the `main.cpp` file.

## License
This is an open source project under the GNU GPLv3.


---

## Usage

### Set global variables in the file `data.h`:

Specify:
 - the number of variables in `const unsigned int n`. 
If the goal is to fit a dataset and generate data from it, then this number can be smaller or equal to the number of variables in the input dataset. If this number is smaller, then only the first `n` bits (per line) of the dataset will be read (starting from the left);
 - the location of the input directory in `const string INPUT_directory`: all the input files must be placed in that folder;
 - the location of the output directory in `const string OUTPUT_directory`: all generated files will be placed in that folder.

### INPUT and OUTPUT files:

Input files must be stored in the INPUT folder. Depending on the functions you are using, you must provide the following input files:
 - **a binary datafile**, if the model is fitted from a given dataset; The input file should be written as strings of `0`'s and `1`'s, with one datapoint per line; see example in `INPUT/Ex_n15_K10_RandM0_dataset_N1000.dat`;
 - **a model file**, if the chosen model is a pairwise model that is not fully connected; Interactions of the fitted models are provided through a file that contains one interaction per line; interactions can be provided in a binary format, with a `1` for each variable included in the interaction and a `0` for the other variable; see example in `INPUT/Ex_PairwiseModel_n4.dat`;
 - **a "matrix file"**, if you want to generate data with a specific values of a) the moments, or b) the covariance, or c) the Pearson correlation coefficients;
    - Specifying the moments (see example in `INPUT/Matrix/Ex_Moments_n4_Bin.dat`):



