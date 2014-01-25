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

/**Eigenvalue solver using LAPACKe
 */
class LapackeEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve, 
								bool assureEigenOrthonormality = true
								);
  static double AssureEigenOrthonormality(EigenInformation * eigenData
										  ); ///Returns maximum deviation in square.

  static void RescaleEigenvectors(EigenInformation * eigenData
									);

private:
  LapackeEigenvalueSolver() {}
};


#endif
