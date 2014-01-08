#include "WorkerData.hh"

WorkerData::WorkerData(CMatrix * _HamiltonianMatrix,
					   ParametrizedCurve * _myCurve,
					   Potential * _myPotential,
					   vector<BasisFunction> _myBasisFunctions,
					   unsigned int _numberOfGLPoints, 
					   double _hbarTimesLambda, 
					   double _massOverLambda2, 
					   double _couplingCoefficient,
					   unsigned int _m1,
					   unsigned int _m2,
					   unsigned int _n1,
					   unsigned int _n2 
					   )
  :HamiltonianMatrix(_HamiltonianMatrix),
   myCurve(_myCurve),
   myPotential(_myPotential),
   myBasisFunctions(_myBasisFunctions),
   numberOfGLPoints(_numberOfGLPoints),
   hbarTimesLambda(_hbarTimesLambda),
   massOverLambda2(_massOverLambda2),
   couplingCoefficient(_couplingCoefficient),
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
