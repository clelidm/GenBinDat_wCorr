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
 - **a "matrix file"**, if you want to generate data with a specific values of a) the moments, or b) the covariance, or c) the Pearson correlation coefficients; see below for details.
 
 #### a) **Specifying the matrix of moments** (see example in `INPUT/Matrix/Ex_Moments_n4_Bin.dat`)
 
The file must contain 1rst and 2nd order moments of the binary variables. The binary variables `Si` must take values `0` or `1`, so that the moments are computed as:
 - 1rst order: `<Si> = P[Si=1]` in the data, which is the probability that `Si` is equal to `1` in the data (for neuronal data, this is similar to the firing rate of the neurons);
 - 2nd order:  `<Si Sj> = P[Si=1 and Sj=1]`,  which is the probability that `Si` and `Sj` are both equal to `1` in the data;

The input file should be written following the format of the example file: `INPUT/Matrix/Ex_Moments_n4_Bin.dat`

*Important:* The value of the moments can only variate between `0` and `1`.
A moment equal to `0.5` corresponds to an unbiased observable (e.g., probability that `Si=1` is `0.5`)
	 whereas a moment equal to `0` or `1` corresponds to an extreme biased observable (e.g., probability that `Si=1` is `0` or `1`),
	 which are impossible to reproduce with the probabilistic model (they will give you parameters with infinite values)
	 ==> don’t use exactly `0` or `1`

#### b) **Specifying the covariance matrix** (see example in `INPUT/Matrix/Ex_Cov_n4_Bin.dat`)

The file must contain 1rst order moments of the binary variables and the coefficients of the covariance matrix.
The binary variables `Si` should take values `0` or `1`, so that these quantities are computed as:
 - 1rst order moment: `<Si> = P[Si=1]` in the data,  i.e. the probability that `Si` is equal to `1` in the data (for neuronal data, this is similar to the firing rate of the neurons);
 - Covariance:  `Cov(i,j) = <Si Sj> - <Si> <Sj>` , where `<Si Sj>` are the 2nd order moments defined as:
	 	 	 	 `<Si Sj> = P[Si=1 and Sj=1]`,  which is the probability that `Si` and `Sj` are both equal to `1` in the data;

The input file should be written following the format of the example file: `INPUT/Matrix/Ex_Cov_n4_Bin.dat`

*Important:* Use realistic values of the covariance. The value of the moments (i.e., `<Si>` and `<Si Sj>`) can variate between `0` and `1`,
You must take covariance values such that this is respected.
A moment equal to `0.5` corresponds to an unbiased observable (e.g., probability that `Si=1` is `0.5`)
	 whereas a moment equal to 0 or 1 corresponding to an extreme biased observable (e.g., probability that `Si=1` is `0` or `1`),
	 which are impossible to reproduce with the probabilistic model (they will give you parameters with infinite values)
 ==> don’t use exactly `0` or `1`
 
#### c) **Specifying the Pearson correlation matrix** (see example in `INPUT/Matrix/Ex_Cov_n4_Bin.dat`)

The file must contain 1rst order moments of the binary variables and the correlation coefficients. The binary variables `Si` should take values `0` or `1`, so that the quantities are computed as:
 - 1rst order moment: `<Si> = P[Si=1]` in the data,  i.e. the probability that `Si` is equal to `1` in the data (for neuronal data, this is similar to the firing rate of the neurons);
 - Correlation Coeff:  `Corr(i,j) = (<Si Sj> - <Si> <Sj>)/sig_i/sig_j` , where `<Si Sj>` are the 2nd order moments defined as:
	 	 	 	 `<Si Sj> = P[Si=1 and Sj=1]`,  which is the probability that `Si` and `Sj` are both equal to `1` in the data;
	 	 and where `sig_i` is the standard deviation of `Si` defined as:
	 	 	 	 `sig_i = std(<Si^2> - <Si>^2)`, which is equal to `sig_i = std(<Si> - <Si>^2)` for the binary variable `Si` in `{0,1}`

The input file should be written following the format of the example file: `INPUT/Matrix/Ex_Cov_n4_Bin.dat`

*Important:* Use realistic values of Pearson correlation. Pearson correlation coefficients can take any values between `-1` and `1`.
The value of the moments (i.e., `<Si>` and `<Si Sj>`) can variate between `0` and `1`,
A moment equal to `0.5` corresponds to an unbiased observable (e.g., probability that `Si=1` is `0.5`)
	 whereas a moment equal to `0` or `1` corresponding to an extreme biased observable (e.g., probability that `Si=1` is `0` or `1`),
	 which are impossible to reproduce with the probabilistic model (they will give you parameters with infinite values)
 ==> don’t use exactly `0` or `1`
Pearson correlation coefficient exactly equal to `-1` or `1` are also not reproducible

## Available functions for generating data with a given matrix of moments, covariance, or correlations

### Read and store the input dataset:
The function map<uint32_t, unsigned int> read_datafile(unsigned int *N) reads the dataset available at the location specified in the variable const string datafilename in data.h. The dataset is then stored in the a structure map<uint32_t, unsigned int> Nset that map each observed states to the number of times they occur in the dataset.

### Fit the data with a model:

### Generate data from the model:

## Available functions for generating data with same firing rate and correlation structure than an input dataset

### Read the matrix of moments, covariance or correlations

### Fit the data with a model:

After fitting, the fitted parameters are available in the list of interactions `list_I` that was given as an argument to the Boltzmann_learning function.
You can print the values of these parameters in the terminal using the function `void PrintTerm_ListInteraction(list<Interaction> list_I)`.

### Generate data from the model:
Once the model is fitted, you can generate data from it using one of the two following functions:
 - The function `void Sample_dataset(list<Interaction> list_I, string output_filename, unsigned int N=1000)` samples a dataset with `N` datapoints from the fitted model defined in `list_I`.
 - The function void Sample_dataset_AND_Print_ModelData_Info(list<Interaction>& list_I, string output_filename, unsigned int N=1000) does the same, while also filling in information about the model and the data in the Interactions of list_I. More precisely, the function computes the model and data averages of the operators of the model and fill this information in list_I (i.e., respectively in I.av_M and I.av_D for each interaction I).
This functions are defined in the file `Generate_data_exact.cpp`.

