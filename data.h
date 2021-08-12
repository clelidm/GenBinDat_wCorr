#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/******************    To Adapt: at each run    *********************/
/********************************************************************/
  const unsigned int n = 4;    // number of spins: only for function "print_basis" to use "bitset<n>()"

  // FOR INPUT FILES:
  const string INPUT_directory = "./INPUT/";

  // FOR OUTPUT FILES:
  const string OUTPUT_directory = "./OUTPUT/";

/********************************************************************/
/**************************    CONSTANTS    *************************/
/*************************    Don't modify   ************************/
/********************************************************************/
  const uint32_t un = 1;
  const uint32_t NOp_tot = (un << n) - 1;     // total number of operators = 2^n - 1

/********************************************************************/
/**************************    STRUCTURES    ************************/
/********************************************************************/
struct Interaction
{
  uint32_t Op;      // binary operator associated to the interaction
  unsigned int k;   // order of the interaction

  double g_Ising;   // parameter of the interaction in {-1,+1} representation  
  double g_Bin;     // parameter of the interaction in {0,1} representation

  double av_D;      // empirical average ("D" for data)
  double av_M;      // average in the Model
};

struct CorrelationM
{
  uint32_t Op;      // binary operator associated to the interaction
  double av_D_Bin;        // empirical average ("D" for data) in the {0, 1}-convention
  double CorrelationCoeff_Bin;  // Correlation Coefficient in the {0, 1}-convention  
};

