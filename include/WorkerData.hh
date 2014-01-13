#ifndef WorkerData_hh
#define WorkerData_hh 1

#include "Globals.hpp"
#include "Matrix.hpp"
#include "ParametrizedCurve.hh"
#include "Potential.hh"
#include <vector>
#include "BasisFunction.hh"
#include "RLException.hh"
#include "HarmonicBasisFunction.hh"


using namespace std;


/**
Class containing the data used by each worker thread.
Essential for the use in RLLib.
*/
class WorkerData
{
public:
  WorkerData(CMatrix * _HamiltonianMatrix, ///Pointer to the Hamilton matrix in use.
			 ParametrizedCurve * _myCurve, ///Pointer to the ParametrizedCurve in use.
			 Potential * _myPotential, ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
			 vector<BasisFunction> _myBasisFunctions, ///Basis functions in use.
			 uint _numberOfGLPoints, ///Number of GL points in use.
			 double _hbarTimesLambda, ///Used for k-E transformation.
			 double _massOverLambda2, ///Used for k-E transformation.
			 double _couplingCoefficient, ///Used only in 2-particle computations.
			 HarmonicBasisFunction * _harmonicBasisFunction,
			 uint _m1, /// Specifies the submatrix to use.
			 uint _m2, /// Specifies the submatrix to use.
			 uint _n1, /// Specifies the submatrix to use.
			 uint _n2 /// Specifies the submatrix to use.
			 ); ///Constructor, basically initialize all values.

  ~WorkerData();
  
  CMatrix * HamiltonianMatrix; ///Pointer to the Hamilton matrix in use.
  ParametrizedCurve * myCurve; ///Pointer to the ParametrizedCurve in use.
  Potential * myPotential; ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
  vector<BasisFunction> myBasisFunctions;
  uint numberOfGLPoints; ///Number of GL points in use.
  double hbarTimesLambda; /// Used for k-E transformation.
  double massOverLambda2; /// Used for k-E transformation.
  double couplingCoefficient; ///The coupling coefficient in 2-particle interactions.
  HarmonicBasisFunction * myHarmonicBasisFunction; ///Harmonic basis function if we are integrating harmonically.
  uint m1; /// Specifies the submatrix to use.
  uint m2; /// Specifies the submatrix to use.
  uint n1; /// Specifies the submatrix to use.
  uint n2; /// Specifies the submatrix to use.
};




#endif
