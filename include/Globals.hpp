#ifndef Globals_hpp
#define Globals_hpp 1


#include "lapacke.h"
#include "lapacke_utils.h"


using namespace std;


#include <complex>
typedef complex<double> ComplexDouble;


#include "Matrix.hpp"
typedef Matrix<ComplexDouble> CMatrix;


#include "SimpsonIntegrator.hpp"
typedef SimpsonIntegrator<ComplexDouble> Integrator;

#define HBAR 1
#define MASS 1

#define foru(var, to)		\
  for(int var = 0; var < to; ++var)

#define foruu(var1, var2, to) \
  for(int var1 = 0; var1 < to; ++var1) \
	for(int var2 = 0; var2 < to; ++var2) \

#define SGN(x) ( (x > 0) - (x < 0) )

#define foric(type, object, it)											\
  for(type::const_iterator it = object.begin(); it!=object.end(); ++it)

//#define make_cplex(re, im) lapack_make_complex_double(re, im)

#define XMIN (-20.)
#define XMAX (20.)
#define NUMBER_OF_X_SAMPLE_POINTS 2000
#define K_CUTOFF 50.
#define NUMBER_OF_K_VALUES 100
#define TRIANGLE_KMID 2.0
#define TRIANGLE_KDEPTH 2.0
  

#endif
