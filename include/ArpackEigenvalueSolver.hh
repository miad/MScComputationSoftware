#ifndef ArpackEigenvalueSolver_hh
#define ArpackEigenvalueSolver_hh 1

#ifdef USE_MKL_LAPACKE
#include "mkl_lapacke.h"
#else
#include "lapacke.h"
#include "lapacke_utils.h"
#endif
#include "ardnsmat.h"
#include "arcomp.h"
#include "ardscomp.h"

#include <vector>
#include <iostream>
#include "Globals.hpp"
#include "Matrix.hpp"
#include "RLException.hh"
#include "CommandLineArgument.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "EigenInformation.hh"


typedef arcomplex<double> ARComplexDouble;
typedef double ARfloat;

using namespace std;

/**Eigenvalue solver using ARPACKpp
 */
class ArpackEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve,
								  uint numberOfEigenvalues,
								  ComplexDouble shift = 0.0,
								  bool findEigenvectors = true
								);
private:
  ArpackEigenvalueSolver() {}
};


#endif
