#ifndef LapackeEigenvalueSolver_hh
#define LapackeEigenvalueSolver_hh 1

//#define LAPACK_BUBBLE_SORT 1

#include "lapacke.h"
#include "lapacke_utils.h"
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "Globals.hpp"
#include "Matrix.hpp"
#include "RLException.hh"
#include "CommandLineArgument.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "EigenInformation.hh"

using namespace std;

/**Eigenvalue solver using LAPACKe
 */
class LapackeEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve,
								  bool ComputeEigenvectors = true
								);
private:
  LapackeEigenvalueSolver() {}
};


#endif
