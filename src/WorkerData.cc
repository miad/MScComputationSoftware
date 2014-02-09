#include "WorkerData.hh"

OneParticleWorkerData::OneParticleWorkerData(CMatrix * _HamiltonianMatrix,
											 ParametrizedCurve * _myCurve,
											 Potential * _myPotential,
											 vector<BasisFunction> _myBasisFunctions,
											 uint _numberOfGLPoints, 
											 double _hbarTimesLambda, 
											 double _massOverLambda2, 
											 HarmonicBasisFunction * _myHarmonicBasisFunction,
											 uint _m1,
											 uint _m2,
											 uint _n1,
											 uint _n2 
											 )
  :HamiltonianMatrix(_HamiltonianMatrix),
   myCurve(_myCurve),
   myPotential(_myPotential),
   myBasisFunctions(_myBasisFunctions),
   numberOfGLPoints(_numberOfGLPoints),
   hbarTimesLambda(_hbarTimesLambda),
   massOverLambda2(_massOverLambda2),
   myHarmonicBasisFunction(_myHarmonicBasisFunction),
   m1(_m1),m2(_m2),n1(_n1),n2(_n2)
{ 
  // Some basic checks.
  if(m1 > m2)
	throw RLException("Tried to initialize OneParticleWorkerData with m1 > m2");
  if(n1 > n2)
	throw RLException("Tried to initialize OneParticleWorkerData with n1 > n2");
  if(myBasisFunctions.empty())
	throw RLException("Tried to initialize OneParticleWorkerData with empty basis function list.");
  if(myCurve == NULL || myPotential == NULL || HamiltonianMatrix == NULL)
	{
	  throw RLException("Tried to create worker data object with a pointer set to NULL.");
	}

}

OneParticleWorkerData::~OneParticleWorkerData()
{
  
}



TwoParticleWorkerData::TwoParticleWorkerData(CMatrix * _HamiltonianMatrix,
											 const PrecomputedInteractionEvaluator * _myPrecomputedInteractionEvaluator,
											 uint _m1,
											 uint _m2,
											 uint _n1,
											 uint _n2 
											 )
  :HamiltonianMatrix(_HamiltonianMatrix),
   myPrecomputedInteractionEvaluator(_myPrecomputedInteractionEvaluator),
   m1(_m1),m2(_m2),n1(_n1),n2(_n2)
{ 
  // Some basic checks.
  if(m1 > m2)
	throw RLException("Tried to initialize TwoParticleWorkerData with m1 > m2");
  if(n1 > n2)
	throw RLException("Tried to initialize TwoParticleWorkerData with n1 > n2");
  if(myPrecomputedInteractionEvaluator == NULL)
	{
	  throw RLException("Tried to create worker data object with PrecomputedInteractionEvaluator pointer set to NULL.");
	}
  if(HamiltonianMatrix == NULL)
	{
	  throw RLException("Tried to create worker data object with HamiltonianMatrix pointer set to NULL.");

	}
}

TwoParticleWorkerData::~TwoParticleWorkerData()
{
  
}
