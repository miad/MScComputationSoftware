#ifndef LapackeEigenvalueSolver_hh
#define LapackeEigenvalueSolver_hh 1

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



class LapackeEigenvalueSolver
{
public:
  static EigenInformation Solve(CMatrix * toSolve);
private:
  LapackeEigenvalueSolver() {}
};


#endif
