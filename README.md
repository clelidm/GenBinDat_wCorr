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

### Read the matrix of moments, covariance or correlations

To create a spin model that can reproduce specific values of first and second moments of the binary variables, we first generated a fully connected pairwise model and then fits its parameters to reproduce the specified moments.

### Fit the data with a model:

To fit a fully connected pairwise model defined in `list_I` on the chosen values of the 1rst and 2nd order moments (same for covariance or correlations), use the function `double BoltzmannLearning_Bin(list<Interaction> &list_I, unsigned int N)`. The argument `unsigned int N` is the number of datapoints in the dataset (real or fictive) that was used to generated the moments/covariance matrix/correlations.

After fitting, the fitted parameters are stored in the list of interactions `list_I` that was given as an argument to the Boltzmann learning function.
You can print the values of these parameters in the terminal using the function `void PrintTerm_ListInteraction(list<Interaction> list_I)`.

### Generate data from the model:
Once the model stored in `list_I` is fitted, there are two functions available to generate data:
 - the function `void SampleData_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1000)` samples a dataset with `N` datapoints from the fitted model defined in `list_I` and print the dataset in `outputfilename`;
 - the function `void SampleData_AND_SaveDataInfo_Bin(list<Interaction>& list_I, string outputfilename, unsigned int N=1000)` does the same, while also filling in information about the model and the generated data in `list_I`. More precisely, the function computes the model and data averages of the operators of the model and fill this information in `list_I` (i.e., respectively in `I.av_M` and `I.av_D` for each interaction `I`).
You can then print this information in the terminal using `void PrintTerm_ModelDataInfo_Bin(list<Interaction> list_I, unsigned int N=1)`, or in an outputfile specified in `outputfilename` using `void PrintFile_ModelDataInfo_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1)`, where `N` is the number of datapoints of the generated dataset.

This functions are defined in the file `Generate_data_exact.cpp`.

## Available functions for generating data with same firing rate and correlation structure than an input dataset

### Read and store the input dataset:
The function `map<uint32_t, unsigned int> read_datafile(string datafilename, unsigned int *N)` reads the dataset provided as an argument in the variable  `string datafilename`. See section `INPUT and OUTPUT file` of this document for information on the format of the input datafile.
The argument `unsigned int *N` is a pointer in which the function will store the total number of datapoints found in the datafile (see example of usage in `main.cpp`). 

The dataset is then stored in the a structure `map<uint32_t, unsigned int> Nset` that map each observed states to the number of times they occur in the dataset. 
You can print this information using the function `void PrintFile_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename)`, where you can specify the name of the output file as an argument in `string OUTPUTfilename`.

### Define a spin model:

A spin model is stored in a list of `Interaction`:  `list<Interaction>`.
For more information, the structure `Interaction` is defined in `data.h`. 
Each `Interaction I` contains the following attributes:
 - the spin operator associated to that interaction, stored as an integer in `I.Op`;
 - the value of the real parameter `I.g` associated to the interaction, stored in `double I.g`;
 - the value of the model average of the operator `I.Op`, stored in `I.av_M`;  --> this value is initially set to `0`, and can be computed after definition of the model 
 - the value of the empirical average of the operator `I.Op`, stored in `I.av_D` --> this value is initially set to `0`, and can be computed when a dataset is generated.

Throughout the program, models (i.e., the list of operators and parameters of the model) are often named starting with `list_I`.

The program offers different ways to define a spin model (i.e., a list of interactions `list<Interaction>`).  These functions are defined in the file `Models_Ex.cpp`. Here is a list:

#### Pre-defined spin models:
 - **Independent model**: the function `list<Interaction> IndepModel(map<uint32_t, unsigned int> Nset, unsigned int N, double *LogLi)` creates an independent model with one field on each spin variable; the parameters of the models are directly fitted to the dataset with `N` datapoints stored in `Nset`.
 - **Fully connected pairwise model**: the function `list<Interaction> FullyConnectedPairwise()` creates a fully connected pairwise model, i.e., with a field on each spin and all the pairwise interactions; the resulting model has `K=n(n+1)/2` interactions; Note that the value of the parameters are not yet specified; this model must be fitted to a dataset before being used to generated data (see section "Fit the data with a model" below).

#### Pairwise spin models specified by the user through an input file: (see example 1.c. in the `int main()` function)
Specific models can be uploaded through an input file. The input file must have the following form:
 - the operator of the model must be written in the first column, in one of these two formats:
   - (a) as a binary representation of the spin involved in the interaction; 
   - (b) as the integer value of that binary representation;
 - operators can only be of order 1 or 2 (i.e., one-body or two-body interactions).

Here are some examples of the two representations of an operator in a 4-spin system:  
 - a field operator on the last digit would be encoded as the binary number `0001`, which has the decimal integer value `1`  -->   0001 = 1
 - a pairwise operator s1 and s2 would be written as the binary number `0011`,  which has the integer representation `3`  -->   0011 = 3

See ``INPUT/Ex_PairwiseModel_n4.dat`` for an example of a file written in format (a).

Two read these files from the `int main()` function, use the function:
 - `list<Interaction> Read_Op_BinaryRepresentation_fromFile(string datafilename);` if the operators are written under binary format (a);
 - `list<Interaction> Read_Op_IntegerRepresentation_fromFile(string datafilename);` if the operators are written under integer format (b);

### Fit the data with a model:
To fit a model previously defined in `list_I` on the dataset stored in `Nset`, use the function `double BoltzmannLearning_Bin_FromData(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N)`. The argument `unsigned int N` is the number of datapoints in the original datafile, that is computed by the function `read_datafile` (see example in `main.cpp`).

After fitting, the fitted parameters are stored in the list of interactions `list_I` that was given as an argument to the Boltzmann learning function.
You can print the values of these parameters in the terminal using the function `void PrintTerm_ListInteraction(list<Interaction> list_I)`.

### Generate data from the model:
Once the model stored in `list_I` is fitted, there are two functions available to generate data:
 - the function `void SampleData_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1000)` samples a dataset with `N` datapoints from the fitted model defined in `list_I` and print the dataset in `outputfilename`;
 - the function `void SampleData_AND_SaveDataInfo_Bin(list<Interaction>& list_I, string outputfilename, unsigned int N=1000)` does the same, while also filling in information about the model and the generated data in `list_I`. More precisely, the function computes the model and data averages of the operators of the model and fill this information in `list_I` (i.e., respectively in `I.av_M` and `I.av_D` for each interaction `I`).
You can then print this information in the terminal using `void PrintTerm_ModelDataInfo_Bin(list<Interaction> list_I, unsigned int N=1)`, or in an outputfile specified in `outputfilename` using `void PrintFile_ModelDataInfo_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1)`, where `N` is the number of datapoints of the generated dataset.

This functions are defined in the file `Generate_data_exact.cpp`.
