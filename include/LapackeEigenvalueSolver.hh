#ifndef LapackeEigenvalueSolver_hh
#define LapackeEigenvalueSolver_hh 1

#include <iostream>
#include "MKL_options.hpp"
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
