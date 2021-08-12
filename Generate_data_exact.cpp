#include <list>
#include <fstream>
#include <map>
#include <bitset> 

using namespace std;
/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/*********************   OPERATOR EVALUATED on a STATE   **********************/
/***********************   depend on the convention    ************************/
/******************************************************************************/

/****************************   !! Convention !!   ****************************/
/****     For Ising:      {1 in data <--> -1 in model}                     ****/
/****                     {0 in data <--> +1 in model}                     ****/
/******************************************************************************/
/****     For Bin:        {1 in data <--> 1 in model}                      ****/
/****                     {0 in data <--> 0 in model}                      ****/
/******************************************************************************/

/******************************************************************************/
/************* PARTITION FUNCTION and PROBA of ALL STATES *********************/
/******************************************************************************/
// return the probability (not normalised) of all the states and compute the partition function
double* Probability_AllStates_Bin(list<Interaction> list_I, double *Z);    // Convention s_i in {0;1} 

/******************************************************************************/
/**************************     SAMPLE DATASET    *****************************/
/******************************************************************************/
void SampleData_Aux(double *p, double Z, string outputfilename, unsigned int N)
{
// Compute Cumulative distribution:
  double *cumul = (double*)malloc((NOp_tot+1)*sizeof(double));

  cumul[0]=p[0]/Z;
  //cout << 0 << "\t " << cumul[0] << endl;

  for(unsigned int i=1; i <= NOp_tot; i++)
  {
    cumul[i] = cumul[i-1] + p[i]/Z;
    //cout << i << "\t " << cumul[i] << endl;
  }
  //cout << cumul[NOp_tot] << endl;

// Randomly sampled "N" data points using the Cumulative:
  double eps=0;
  uint32_t j=0;

// OUTPUT FILE:
  fstream file((outputfilename).c_str(), ios::out);

  for(int i=0; i<N; i++)
  {
    eps=drand48();
    j=0;
    while(cumul[j]<eps)
    { j++;  }
    file << bitset<n>(j) << endl;
  }

  file.close();
  free(cumul);
}

void SampleData_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1000)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *P=Probability_AllStates_Bin(list_I, &Z);

  SampleData_Aux(P, Z, outputfilename, N);
  free(P);
}

/******************************************************************************/
/******************************************************************************/
/**************************     GET INFO      *********************************/
/********************    ON THE GENERATED DATA    *****************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/**************************     MODEL AVERAGES    *****************************/
/******************************************************************************/
void Model_averages_Bin_Aux(double *P, double Z, list<Interaction> &list_I);

void Model_averages_Bin(list<Interaction> &list_I)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *P=Probability_AllStates_Bin(list_I, &Z);
  Model_averages_Bin_Aux(P, Z, list_I);
  free(P);
}

/******************************************************************************/
/*************************     EMPIRICAL AVERAGES     *************************/
/******************************************************************************/
void empirical_averages_Bin(map<uint32_t, unsigned int> Nset, list<Interaction> &list_I, unsigned int N);

/******************************************************************************/
/**************************     SAMPLE DATASET    *****************************/
/*************  and fill info on model + dataset    ***************************/
/********  i.e. compute Model averages and data averages    *******************/
/***************  and print them with the model information    ****************/
/******************************************************************************/
map<uint32_t, unsigned int> SampleData_AND_SaveDataInfo_Aux(double *P, double Z, string outputfilename, unsigned int N)
{
// Compute Cumulative distribution:
  double *cumul = (double*)malloc((NOp_tot+1)*sizeof(double));

  cumul[0]=P[0]/Z;

  for(int i=1; i <= NOp_tot; i++)
  {
    cumul[i] = cumul[i-1] + P[i]/Z;
  }

// Randomly sampled "N" data points using the Cumulative:
  double eps=0;
  uint32_t datapt=0;

// OUTPUT FILE:
  fstream file((outputfilename).c_str(), ios::out);

// ***** data is also stored in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set

  for(unsigned int i=0; i<N; i++)
  {
    eps=drand48();
    datapt=0;

    while(cumul[datapt]<eps)
    { datapt++;  }

    file << bitset<n>(datapt) << endl;
    Nset[datapt] += 1;
  }

  file.close();
  free(cumul);

  return Nset;
}

void SampleData_AND_SaveDataInfo_Bin(list<Interaction>& list_I, string outputfilename, unsigned int N=1000)
{
// Compute un-normalized probability distribution and Partition function Z:
  double Z=0;
  double *P=Probability_AllStates_Bin(list_I, &Z);

  map<uint32_t, unsigned int> Nset = SampleData_AND_SaveDataInfo_Aux(P, Z, outputfilename, N);

  // Model averages:
  Model_averages_Bin_Aux(P, Z, list_I);
  free(P);

  // Empirical averages (from the generated dataset):
  empirical_averages_Bin(Nset, list_I, N);
}

/******************************************************************************/
/**************************   PRINT MODEL IN TERMINAL   ***********************/
/*****************************    ALL INTERACTIONS    *************************/
/******************************************************************************/
double std_M(double av_M, double Nd)
{
  return (Nd==1) ? 0 : sqrt( (1.+av_M)*(1-av_M)/Nd )/2.;
}

void PrintTerm_ModelDataInfo_Bin(list<Interaction> list_I, unsigned int N=1)
{
  list<Interaction>::iterator it;
  double Nd = (double) N;

  cout << "# In case a dataset was generated: N=" << N << endl;
  cout << "# 1:Op \t 2:order \t 3:g_Ising \t 4:Model_av \t 5:Data_av \t 6:Std_of_avD_from_avM(for dataset with N datapoints)" << endl;
  for (it=list_I.begin(); it!=list_I.end(); it++)
  {
    cout << "Op = " << (*it).Op << " = " << bitset<n>((*it).Op)  << "\t " << bitset<n>((*it).Op).count() << "\t g_Bin = " << (*it).g_Bin;
    cout << "\t Model_av = " << (*it).av_M;
    cout << "\t Data_av = " << (*it).av_D;
    cout << "\t Std_of_avD_from_avM = " << std_M((*it).av_M, Nd) << endl;
  }
  cout << endl;
}

void PrintFile_ModelDataInfo_Bin(list<Interaction> list_I, string outputfilename, unsigned int N=1)
{
  list<Interaction>::iterator it;
  double Nd = (double) N;

  fstream file((outputfilename).c_str(), ios::out);

  file << "# In case a dataset was generated: N=" << N << endl;
  file << "# 1:Operator \t 2:Order \t 3:g_Bin \t 4:Model_average \t 5:Data_average";
  file << "\t 6:Std_of_avD_from_avM(for dataset with N datapoints)";
  file << "\t 7:Operator_bit_representation" << endl;
  file << "### " << endl;
  for (it=list_I.begin(); it!=list_I.end(); it++)
  {
    file << (*it).Op << "\t " << bitset<n>((*it).Op).count() << "\t " << (*it).g_Bin;
    file << "\t " << (*it).av_M << "\t " << (*it).av_D;
    file << "\t " << std_M((*it).av_M, Nd);
    file << "\t " << bitset<n>((*it).Op) << endl;
  }
  file << endl;

  file.close();
}

