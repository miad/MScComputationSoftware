#ifndef CompositeBasisFunction_hh
#define CompositeBasisFunction_hh 1

#include "Globals.hpp"
#include "RLException.hh"
#include <stdio.h>
#include <ctype.h>
#include <cstring>
#include <map>
#include <iostream>
#include <algorithm>
#include "HermiteEvaluator.hh"
#include "BasisFunction.hh"
#include "SpecificUnits.hh"
#include "EigenInformation.hh"
#include "HarmonicBasisFunction.hh"
#include "ParametrizedCurve.hh"

using namespace std;

/**Used in 2-particle systems for interaction.
 */
class CompositeBasisFunction
{
public:
  CompositeBasisFunction(vector<BasisFunction> _functions, ///Basis functions to use.
						 EigenInformation * _myInformation, ///Eigeninformation. Ownership is not taken. May be deleted after constructor call.
						 const ParametrizedCurve * _KCurve
						 ); ///Constructor. Ownership of the objects are NOT passed, and they are NOT used after the constructor.

  CompositeBasisFunction(const HarmonicBasisFunction * _myHarmonicBasisFunction,
						 EigenInformation * myInformation
						 ); ///For using with HarmonicOverride...

  ~CompositeBasisFunction(); ///Destructor.

  ComplexDouble Eval(const double & x, ///x-value for evaluation.
					 uint pIndex ///Basis function index.
					 ); ///Evaluate the function using parameter x.

  ComplexDouble GetE(uint pIndex ///The index of evaluation.
					 ) const; ///Return the energy corresponding to p-index.

  uint GetSize() const; ///Number of eigenvalues / eigenvectors.

  bool IsHarmonic() const; ///Check if the object was initialized as harmonic object or not.

protected:
  double HarmonicEval(const double & x, ///x-value.
					  uint pIndex ///Basis function index.
					  ) const; ///Evaluate.

  ComplexDouble KEval(const double & x, ///x-value
					  uint pIndex ///Basis function index
					  );  ///Evaluate using parameter x.

private:
  vector<BasisFunction> functions; ///Basis functions, in the event of non-harmonic. Otherwise empty.
  const HarmonicBasisFunction * myHarmonicBasisFunction; ///Basis functions, in the event of harmonic. Otherwise NULL.

  const EigenInformation * myInformation; ///Eigen information. 
  const ParametrizedCurve * KCurve; ///K-curve, in the event of non-harmonic. Otherwise NULL.
};

#endif
