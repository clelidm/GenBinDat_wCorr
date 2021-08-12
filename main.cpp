//g++ -I /usr/local/include/eigen3/ -I ./LBFGSpp-master/include -O2 main.cpp Models.cpp ReadDataFile.cpp BoltzmannLearning.cpp FitFromCorrMatrix.cpp Generate_data_exact.cpp
//time ./a.out
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <map>
#include <set> 
#include <ctime> // for chrono
#include <ratio> // for chrono
#include <chrono> // for chrono
#include <random>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
map<uint32_t, unsigned int> read_datafile(string datafilename, unsigned int *N);   // O(N)  //N = data set size  //READ DATA and STORE them in Nset
void read_Nset (map<uint32_t, unsigned int> Nset, unsigned int N, string OUTPUTfilename); // PRINT Nset in a file  

/******************************************************************************/
/*************************  MODEL: LIST of INTERACTIONS ***********************/
/***************************   in Models.cpp     ******************************/
/******************************************************************************/

/************************** INDEPENDENT MODEL *********************************/
list<Interaction> IndepModel(map<uint32_t, unsigned int> Nset, unsigned int N, double *LogLi);

/**************************** OTHER MODELS ************************************/
list<Interaction> FullyConnectedPairwise();

/**************************** READ MODEL FROM FILE ****************************/
list<Interaction> Read_Op_IntegerRepresentation_fromFile(string datafilename);
list<Interaction> Read_Op_BinaryRepresentation_fromFile(string datafilename);

/************************* PRINT MODEL in TERMINAL ****************************/
void PrintTerm_ListInteraction(list<Interaction> list_I);

/******************************************************************************/
/**************************  CORRELATION MATRIX  ******************************/
/*************************    computed on data    *****************************/
/******************************************************************************/
list<CorrelationM> CorrelationMatrix(map<uint32_t, unsigned int> Nset, unsigned int N);
void Print_Term_CorrelationMatrix(list<CorrelationM> list_Corr);

/******************************************************************************/
/*************************  MODEL FROM CORRELATION MATRIX *********************/
/**************************   in FitFromCorrMatrix.cpp     ********************/
/******************************************************************************/
//double BoltzmannLearning_Ising_CorrMatrix(string CorrMatrix_File, unsigned int N);
list<Interaction> ModelFrom_MomentMatrix(string MomentMatrix_File, bool *error);
list<Interaction> ModelFrom_CorrMatrix(string CorrMatrix_File, bool *error);

/******************************************************************************/
/**************************** Boltzmann Learning ******************************/
/******************************************************************************/
double BoltzmannLearning_Bin_FromData(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);
double BoltzmannLearning_Bin(list<Interaction> &list_I, unsigned int N);

/******************************************************************************/
/**********************  Generate Data from fitted model **********************/
/******************************************************************************/
// RANDOM GENERATOR
void Initialise_RandGenerator()
{
  int seed = (unsigned)time(NULL);
  srand48(seed);      //for drand48
}

// If parameters are specified for the {0,1}-convention:
// !!! This version only works for models with field and/or pairwise interactions !!
void SampleData_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1000);
void SampleData_AND_SaveDataInfo_Bin(list<Interaction>& list_I, string outputfilename, unsigned int N=1000);

// Print Info about the model and generated dataset:
void PrintTerm_ModelDataInfo_Bin(list<Interaction> list_I, unsigned int N=1);
void PrintFile_ModelDataInfo_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1);

/******************************************************************************/
/************************** MAIN **********************************************/
/******************************************************************************/
//Rem : ---> 2^30 ~ 10^9
int main()
{
//*******************************************
//********** READ DATA FILE:  ***************     -->  data are put in Nset:
//*******************************************    // Nset[mu] = # of time state mu appears in the data set

  cout << "--->> Create OUTPUT Folder: (if needed) ";
  system( ("mkdir -p " + OUTPUT_directory).c_str() );
  cout << endl;
 
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*************************  INFERRING A MODEL FROM DATA:  **********************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "************************************  Read dataset:  **************************************" << endl;

  unsigned int N = 0;  // set in 'read_datafile()'
  string INPUT_Data_Filename = INPUT_directory + "Ex_n15_K10_RandM0_dataset_N1000.dat";
  string OUTPUT_Data_Filename = OUTPUT_directory + "Ex_n15_K10_RandM0_dataset_N1000.dat";

  map<uint32_t, unsigned int> Nset = read_datafile(INPUT_Data_Filename, &N);  //  double Nd = (double) N;
  read_Nset(Nset, N, OUTPUT_Data_Filename + "_Nset.dat");   // print Nset in a file

  cout << endl << "*****************************  Ex 1.a. Independent model:  *********************************" << endl;
  cout << "--->> Inferred field for the Independent Model associated to the dataset '" << INPUT_Data_Filename << "'" << endl << endl;

  double L_indep;  // LogLikelihood of the independent model
  list<Interaction> list_I_indep = IndepModel(Nset, N, &L_indep);
  int K=0; double L=0;

  K = list_I_indep.size();
  L = L_indep;     
  
  cout << "Number of parameters, K = " << K << endl;
  cout << "Max Log-Likelihood, L = " << L << endl;
  cout << "BIC: " << L - K * log( ((double) N) /(2.*M_PI)) / 2. <<  endl << endl;
  PrintTerm_ListInteraction (list_I_indep);

  cout << endl << "*********************  Ex 1.b. Fully connected pairwise model:  *****************************" << endl;

  list<Interaction> list_I_pairwise = FullyConnectedPairwise();
  K = list_I_pairwise.size();

   // ********* TIME TEST with chrono:
  chrono::high_resolution_clock::time_point t1, t2;

  // ************** fit with {0,1}-model:
  t1 = chrono::high_resolution_clock::now(); 
  cout << "Number of parameters to fit: K = " << K << endl;
  cout << "--->> Fitting (gradient descent) for {0,1}-model.." << endl;
  L = BoltzmannLearning_Bin_FromData (Nset, list_I_pairwise, N);

  t2 = chrono::high_resolution_clock::now();
  cout << "Timings test: time machine = " << (chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())/1e9 << " sec." << endl << endl;

  cout << "Max Log-Likelihood, L = " << L << endl;
  cout << "BIC: " << L - K * log( ((double) N) /(2.*M_PI)) / 2. <<  endl << endl;
  PrintTerm_ListInteraction (list_I_pairwise);

  cout << endl << "**************  Ex 1.c. Fit a pairwise model defined from an input file:  *******************" << endl;
  cout << "\t !!  This version only works for models with one-body or two-body interactions (i.e., fields or pairwise interactions) !!" << endl;
// ************** Fit a model defined in input file:
  cout << "Read model from file:" << endl;
  list<Interaction> list_I_Model = Read_Op_BinaryRepresentation_fromFile("INPUT/Ex_PairwiseModel_n4.dat");

  cout << "Number of parameters to fit: K = " << list_I_Model.size() << endl;
  cout << "--->> Fitting (gradient descent) for {0,1}-model.." << endl;
  L = BoltzmannLearning_Bin_FromData (Nset, list_I_Model, N);

  cout << "Max Log-Likelihood, L = " << L << endl;
  cout << "BIC: " << L - K * log( ((double) N) /(2.*M_PI)) / 2. <<  endl << endl;
  PrintTerm_ListInteraction (list_I_Model);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "**************  Ex 2. COMPUTING THE CORRELATION MATRIX of an INPUT DATA FILE:  ************";
  cout << endl << "*******************************************************************************************" << endl;
  // Compute the correlation Coefficient of the dataset in 'INPUT_filename':
  cout << endl << "Correlation matrix of the datafile: " << INPUT_Data_Filename << endl;
  list<CorrelationM> list_Corr = CorrelationMatrix(Nset, N);
  Print_Term_CorrelationMatrix(list_Corr);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************  Ex 3. PAIRWISE MODEL INFERRED FROM A CORRELATION MATRIX:  ***************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "*********************  Ex 3.a. from the matrix of moments:  *******************************" << endl;

  string InputFile_Ex_Moments = "INPUT/Matrix/Ex_Moments_n4_Bin.dat";
  unsigned int N_Matrix = 1000;  // number of datapoints in the datafile for which the correlation matrix is provided 

  cout << endl << "The file should contain 1rst and 2nd order moments of the binary variables." << endl;
  cout << "The binary variables Si should take values `0` or `1`, so that the moments are computed as:" << endl;
  cout << "\t \t 1rst order: <Si> = P[Si=1] in the data,  i.e. the probability that Si is equal to 1 in the data;" << endl;
  cout << "\t \t \t \t \t \t  --> for neuronal data, that would similar to the firing rate;" << endl;
  cout << "\t \t 2nd order:  <Si Sj> = P[Si=1 and Sj=1],  i.e. the probability that Si and Sj are both equal to 1 in the data;" << endl;
  cout << endl << "The input file should be written following the format of the example file: " << InputFile_Ex_Moments << endl;

  cout << "***********" << endl;
  bool error = false;
  list<Interaction> list_I_MomentMatrix_Bin = ModelFrom_MomentMatrix(InputFile_Ex_Moments, &error);
  if (!error)
  {
    L = BoltzmannLearning_Bin(list_I_MomentMatrix_Bin, N_Matrix);
    cout << "Max Log-Likelihood, L = " << L << endl;
    cout << "BIC: " << L - K * log( ((double) N_Matrix) /(2.*M_PI)) / 2. <<  endl << endl;
    PrintTerm_ListInteraction (list_I_MomentMatrix_Bin);
  }

  cout << endl << "*********************  Ex 3.b. from the matrix of Correlations:  ****************************" << endl;

  string InputFile_Ex_Corr = "INPUT/Matrix/Ex_CorrCoeff_n4_Bin.dat";

  cout << endl << "The file should contain 1rst order moments of the binary variables and the correlation coefficients." << endl;
  cout << "The binary variables Si should take values `0` or `1`, so that the quantities are computed as:" << endl;
  cout << "\t \t 1rst order moment: <Si> = P[Si=1] in the data,  i.e. the probability that Si is equal to 1 in the data;" << endl;
  cout << "\t \t \t \t \t \t  --> for neuronal data, that would similar to the firing rate;" << endl;
  cout << "\t \t Correlation Coefficients:  Corr(i,j) = <Si Sj> - <Si> <Sj> , where <Si Sj> are the 2nd order moments defined as:" << endl;
  cout << "\t \t \t \t \t \t <Si Sj> = P[Si=1 and Sj=1],  which is the probability that Si and Sj are both equal to 1 in the data;" << endl;

  cout << endl << "The input file should be written following the format of the example file: " << InputFile_Ex_Corr << endl;

  cout << "***********" << endl;
  list<Interaction> list_I_CorrMatrix_Bin = ModelFrom_CorrMatrix(InputFile_Ex_Corr, &error);
  if (!error)
  {
    L = BoltzmannLearning_Bin(list_I_CorrMatrix_Bin, N_Matrix);
    cout << "Max Log-Likelihood, L = " << L << endl;
    cout << "BIC: " << L - K * log( ((double) N_Matrix) /(2.*M_PI)) / 2. <<  endl << endl;
    PrintTerm_ListInteraction (list_I_CorrMatrix_Bin);
  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "***********************  Ex 4. GENERATING DATA FROM A CHOSEN MODEL:  **********************";
  cout << endl << "*******************************************************************************************" << endl;

  Initialise_RandGenerator();

  cout << endl << "***********  Ex 4. from the model inferred from the matrix of Correlations:  **************" << endl;

  string OutputFile_Ex_DataFromCorr = OUTPUT_directory + "FromCorrelations";
  string OUTPUT_Data = OutputFile_Ex_DataFromCorr + "_Data.dat";
  string OUTPUT_Data_Info = OutputFile_Ex_DataFromCorr + "_Model_AND_Data_Info.dat";

  unsigned int N_new = 1e5;

  cout << "Generate a dataset with N = " << N_new << " datapoints, stored in " << OUTPUT_Data << endl;
  SampleData_AND_SaveDataInfo_Bin(list_I_CorrMatrix_Bin, OUTPUT_Data, N_new);

  cout << "Info about the model and the generated data:" << endl;
  PrintFile_ModelDataInfo_Bin(list_I_CorrMatrix_Bin, OUTPUT_Data_Info, N_new);
  PrintTerm_ModelDataInfo_Bin(list_I_CorrMatrix_Bin, N_new);

  cout << "Correlation matrix of the generated data:" << endl;
  map<uint32_t, unsigned int> Nset_Bin = read_datafile(OUTPUT_Data, &N_new);  //  double Nd = (double) N;
  Print_Term_CorrelationMatrix(CorrelationMatrix(Nset_Bin, N_new));

  OUTPUT_Data = OutputFile_Ex_DataFromCorr + "_Data2.dat";
  cout << "Generate a second dataset with N = " << N_new << " datapoints based on the same correlation matrix." << endl;
  cout << "Dataset is stored in " << OUTPUT_Data << endl << endl;
  SampleData_Bin(list_I_CorrMatrix_Bin, OUTPUT_Data, N_new);

  return 0;
}


