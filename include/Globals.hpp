#ifndef Globals_hpp
#define Globals_hpp 1


#include "lapacke.h"
#include "lapacke_utils.h"


using namespace std;


#include <complex>
typedef complex<double> ComplexDouble;


#include "Matrix.hpp"
typedef Matrix<ComplexDouble> CMatrix;

#include "RLMacros.hpp"


#define HBARC (197.326971812) ///hbar * c in MeV * fm
#define MASSOVERC2 938 ///mass in MeV/c^2

/*#define XMIN (-20.)
#define XMAX (20.)
#define NUMBER_OF_X_SAMPLE_POINTS 2000
#define K_CUTOFF 50.
#define NUMBER_OF_K_VALUES 100
#define TRIANGLE_KMID 2.0
#define TRIANGLE_KDEPTH 2.0
*/

#endif
