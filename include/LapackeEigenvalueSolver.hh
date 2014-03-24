#ifndef LapackeEigenvalueSolver_hh
#define LapackeEigenvalueSolver_hh 1

//#define LAPACK_BUBBLE_SORT 1

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

#ifdef USE_MKL_LAPACKE
#include "mkl_lapacke.h"
#define RIKARD_COMPLEX_TYPE MKL_Complex16
#define RIKARD_COMPLEX_COPY_TO(from, to){				\
	to.real = from.real();								\
	to.imag = from.imag();								\
  }
#define RIKARD_COMPLEX_COPY_FROM(from, to){				\
	to.real() = from.real;								\
	to.imag() = from.imag;								\
  }
#else
#define RIKARD_COMPLEX_TYPE ComplexDouble
#include "lapacke.h"
#include "lapacke_utils.h"
#define RIKARD_COMPLEX_COPY_TO(from, to){		\
	to = from;									\
  }
#define RIKARD_COMPLEX_COPY_FROM(from, to){		\
	to = from;									\
  }
#endif

using namespace std;

/**Eigenvalue solver using LAPACKe
 */
class LapackeEigenvalueSolver
{
public:
  static EigenInformation * Solve(CMatrix * toSolve,
								  bool ComputeEigenvectors = true
								);
  static RIKARD_COMPLEX_TYPE * GetArray(CMatrix * matrix
								  );

private:
  LapackeEigenvalueSolver() {}
};


#endif
