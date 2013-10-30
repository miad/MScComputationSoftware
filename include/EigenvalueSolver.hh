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

using namespace std;


struct EigenInformation
{
  vector<vector<ComplexDouble> > Eigenvectors;
  vector<ComplexDouble> Eigenvalues;
};



class EigenvalueSolver
{
public:
  static EigenInformation Solve(CMatrix * toSolve);
private:
  EigenvalueSolver() {}
};


#endif
