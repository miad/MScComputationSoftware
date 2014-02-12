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
#include "CompositeBasisFunction.hh"
#include "InteractionProperties.hh"
#include "PrecomputedInteractionEvaluator.hh"


using namespace std;


/**
Class containing the data used by each worker thread.
Essential for the use in RLLib.
*/
class OneParticleWorkerData
{
public:
  OneParticleWorkerData(CMatrix * _HamiltonianMatrix, ///Pointer to the Hamilton matrix in use.
			 ParametrizedCurve * _myCurve, ///Pointer to the ParametrizedCurve in use.
			 Potential * _myPotential, ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
			 vector<BasisFunction> _myBasisFunctions, ///Basis functions in use.
			 uint _numberOfGLPoints, ///Number of GL points in use.
			 double _hbarTimesLambda, ///Used for k-E transformation.
			 double _massOverLambda2, ///Used for k-E transformation.
			 HarmonicBasisFunction * _harmonicBasisFunction,
						uint _particleID, ///Particle ID.
			 uint _m1, /// Specifies the submatrix to use.
			 uint _m2, /// Specifies the submatrix to use.
			 uint _n1, /// Specifies the submatrix to use.
			 uint _n2 /// Specifies the submatrix to use.
			 ); ///Constructor, basically initialize all values.

  ~OneParticleWorkerData();
  
  CMatrix * HamiltonianMatrix; ///Pointer to the Hamilton matrix in use.
  ParametrizedCurve * myCurve; ///Pointer to the ParametrizedCurve in use.
  Potential * myPotential; ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
  vector<BasisFunction> myBasisFunctions;
  uint numberOfGLPoints; ///Number of GL points in use.
  double hbarTimesLambda; /// Used for k-E transformation.
  double massOverLambda2; /// Used for k-E transformation.
  HarmonicBasisFunction * myHarmonicBasisFunction; ///Harmonic basis function if we are integrating harmonically.
  uint particleID;
  uint m1; /// Specifies the submatrix to use.
  uint m2; /// Specifies the submatrix to use.
  uint n1; /// Specifies the submatrix to use.
  uint n2; /// Specifies the submatrix to use.
};




/**
Class containing the data used by each worker thread.
Essential for the use in RLLib.
*/
class TwoParticleWorkerData
{
public:
  TwoParticleWorkerData(CMatrix * _HamiltonianMatrix,
						const PrecomputedInteractionEvaluator * _myPrecomputedInteractionEvaluator,
						uint _m1,
						uint _m2,
						uint _n1,
						uint _n2 
						); ///Constructor, basically initialize all values.
  
  ~TwoParticleWorkerData();
  
  CMatrix * HamiltonianMatrix;
  const PrecomputedInteractionEvaluator * myPrecomputedInteractionEvaluator;

  uint m1; /// Specifies the submatrix to use.
  uint m2; /// Specifies the submatrix to use.
  uint n1; /// Specifies the submatrix to use.
  uint n2; /// Specifies the submatrix to use.
};




#endif
