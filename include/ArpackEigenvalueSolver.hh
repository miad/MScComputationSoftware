#ifndef ArpackEigenvalueSolver_hh
#define ArpackEigenvalueSolver_hh 1

#include <stdio.h>
#include <vector>
#include "lapacke.h"
#include "lapacke_utils.h"
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
class ArpackEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve
								);
private:
  ArpackEigenvalueSolver() {}
};


#endif
