#include <fstream>
#include <sstream>
#include <map>
#include <list>
#include <bitset> 

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************   Correlation matrix   ****************************/
/*******************   Compute the Corr Coeffcients   *************************/
/*************************   For a given dataset   ****************************/
/******************************************************************************/
//Number of time an operator is equal to 1 ( = <phi> in the {0,1} representation )
double op_av_Bin(map<uint32_t, unsigned int> Nset, uint32_t op, unsigned int N);

/**********************   COMPUTE ALL THE COEFFICIENTS   **************************/
list<CorrelationM> CorrelationMatrix(map<uint32_t, unsigned int> Nset, unsigned int N)
{
  CorrelationM Corr;
  list<CorrelationM> list_Corr;

  uint32_t Op1 = 1, Op2 = 1, Op = 1;
  double* magne = (double*)malloc(n*sizeof(double));

  // local magnetizations:
  for(int i=0; i<n; i++) // n fields
  {
    Corr.Op = Op1;   
    Corr.av_D_Bin = op_av_Bin(Nset, Op1, N);

    magne[i] = Corr.av_D_Bin;
    Corr.CorrelationCoeff_Bin = Corr.av_D_Bin;

    list_Corr.push_back(Corr);

    Op1 = Op1 << 1;
  }

  // local correlations:
  Op1 = 1;

  for(int i=0; i<n; i++) // n(n-1)/2 pairwise interactions
  {    
    Op2 = Op1 << 1;
    for (int j=i+1; j<n; j++)
    {
      Op = Op1 + Op2;
      
      Corr.Op = Op;
      Corr.av_D_Bin = op_av_Bin(Nset, Op, N);

      Corr.CorrelationCoeff_Bin = Corr.av_D_Bin - magne[i]*magne[j];

      list_Corr.push_back(Corr);

      Op2 = Op2 << 1; 
    }
    Op1 = Op1 << 1;      
  }

  free(magne);

  return list_Corr;
}

/**********************   PRINT MATRIX IN TERMINAL   **************************/
void Print_Term_CorrelationMatrix(list<CorrelationM> list_Corr)
{
  list<CorrelationM>::iterator it_corr;
  cout << "# 1:Op   2:Emp_av_Ising  3:Emp_av_{0,1}   4:CorrCoeff_{0,1}" << endl;

  for(it_corr = list_Corr.begin(); it_corr != list_Corr.end(); it_corr++)
  {
    cout << "Op = " << (*it_corr).Op << " = " << bitset<n>((*it_corr).Op);
    cout << " \t Emp_av_{0,1} = " << (*it_corr).av_D_Bin;

    if(bitset<n>((*it_corr).Op).count() == 2)
    {
      cout << " \t CorrCoeff_{0,1} = " << (*it_corr).CorrelationCoeff_Bin;
    }
    cout << endl;
  }
  cout << endl;
}

/******************************************************************************/
/**************************   FULLY CONNECTED PAIRWISE MODEL   ****************/
/******************************************************************************/
list<double> Read_Input_MatrixFile(string Input_Matrix_File)
{
  list<double> list_MatrixValues;

  ifstream myfile (Input_Matrix_File.c_str());
  string line, temp;    
  stringstream sstream; 
  double av_D = 0;

  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      sstream << line;   // Storing the whole string into string stream
      sstream >> temp;   // extracting first "word" from stream

      if (stringstream(temp) >> av_D)   // Checking the given "word" is double or not
      { 
        //cout << av_D << endl;
        list_MatrixValues.push_back(av_D);  
      }

      temp = "";
      sstream.str("");
    }
    myfile.close();
  }
  //cout << list_MatrixValues.size() << endl;

  return list_MatrixValues;
}

list<Interaction> FullyConnectedPairwise();

list<Interaction> ModelFrom_MomentMatrix(string MomentMatrix_File, bool *error)
{
  list<double> list_MatrixValues = Read_Input_MatrixFile(MomentMatrix_File);
  list<Interaction> list_I = FullyConnectedPairwise();

  //int K = n*(n+1)/2;
  //cout << K << "\t" << list_MatrixValues.size() << "\t" << list_I.size() << endl;

  if (list_MatrixValues.size()!= list_I.size()) 
  { 
    cout << "Error: the number of lines in the Correlation matrix file does not match with ";
    cout << "the number `n` of neurons specified in `data.h`. The number of lines should be n(n+1)/2." << endl; 
    *error=true;
  }
  else
  {
    list<double>::iterator it_av;
    list<Interaction>::iterator it_I;

    it_av = list_MatrixValues.begin();
    for(it_I = list_I.begin(); it_I != list_I.end(); it_I++)
    {
      (*it_I).av_D = (*it_av);
      it_av++;
    }
    *error=false;
  } 

  return list_I;  
}

list<Interaction> ModelFrom_CorrMatrix(string CorrMatrix_File, bool *error)
{
  list<double> list_MatrixValues = Read_Input_MatrixFile(CorrMatrix_File);

  list<Interaction> list_I;

  int K = n*(n+1)/2;
  if (list_MatrixValues.size()!=K) 
  { 
    cout << "Error: the number of lines in the Correlation matrix file does not match with ";
    cout << "the number `n` of neurons specified in `data.h`. The number of lines should be n(n+1)/2." << endl; 
    *error = true;
  }
  else
  {
    *error = false;
    Interaction I;
    I.g_Ising = 0;  I.g_Bin = 0;  I.av_M = 0;

    double* magne = (double*)malloc(n*sizeof(double));
    uint32_t Op1 = 1, Op2 = 1;

    list<double>::iterator it_av = list_MatrixValues.begin();

    I.k = 1;
    for(int i=0; i<n; i++) // n fields
    {
      I.Op = Op1;
      I.av_D = *it_av;
      list_I.push_back(I);    

      magne[i] = I.av_D;

      Op1 = Op1 << 1;
      it_av++;
    }

    Op1 = 1;

    I.k = 2;
    for(int i=0; i<n; i++) // n(n-1)/2 pairwise interactions
    {    
      Op2 = Op1 << 1;
      for (int j=i+1; j<n; j++)
      {
        I.Op = Op1 + Op2; 
        I.av_D = (*it_av) + magne[i]*magne[j];
        list_I.push_back(I);
      
        Op2 = Op2 << 1; 
        it_av++;
      }
      Op1 = Op1 << 1;      
    }
  }

  return list_I;
}
