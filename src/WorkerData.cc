#include "WorkerData.hh"

WorkerData::WorkerData(CMatrix * _HamiltonianMatrix,
					   ParametrizedCurve * _myCurve,
					   Potential * _myPotential,
					   vector<BasisFunction> _myBasisFunctions,
					   uint _numberOfGLPoints, 
					   double _hbarTimesLambda, 
					   double _massOverLambda2, 
					   double _couplingCoefficient,
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
   couplingCoefficient(_couplingCoefficient),
   myHarmonicBasisFunction(_myHarmonicBasisFunction),
   m1(_m1),m2(_m2),n1(_n1),n2(_n2)
{ 
  // Some basic checks.
  if(m1 > m2)
	throw RLException("Tried to initialize WorkerData with m1 > m2");
  if(n1 > n2)
	throw RLException("Tried to initialize WorkerData with n1 > n2");
  if(myBasisFunctions.empty())
	throw RLException("Tried to initialize WorkerData with empty basis function list.");
}

WorkerData::~WorkerData()
{
}
