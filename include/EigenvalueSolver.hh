#ifndef EigenvalueSolver_hh
#define EigenvalueSolver_hh 1

#include "lapacke.h"
#include "lapacke_utils.h"
#include "Globals.hpp"
#include "Matrix.hpp"
#include "RLException.hh"
#include <stdio.h>
#include <vector>
#include "CommandLineArgument.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "EigenInformation.hh"

using namespace std;



class EigenvalueSolver
{
public:
  static EigenInformation Solve(CMatrix * toSolve);
private:
  EigenvalueSolver() {}
};


#endif
