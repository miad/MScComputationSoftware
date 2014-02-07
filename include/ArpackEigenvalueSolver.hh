#ifndef ArpackEigenvalueSolver_hh
#define ArpackEigenvalueSolver_hh 1

#include "lapacke.h"
#include "lapacke_utils.h"
//#include "ardnsmat.h"
//#include "arcomp.h"
//#include "ardscomp.h"

#include <vector>
#include "Globals.hpp"
#include "Matrix.hpp"
#include "RLException.hh"
#include "CommandLineArgument.hh"
#include "CommandLineInterpreter.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "EigenInformation.hh"


//typedef arcomplex<double> ARComplexDouble;
//typedef double ARfloat;

using namespace std;

/**Eigenvalue solver using LAPACKe
 */
class ArpackEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve,
								  uint numberOfEigenvalues,
								  ComplexDouble shift = 0.0
								);
private:
  ArpackEigenvalueSolver() {}
};


#endif
