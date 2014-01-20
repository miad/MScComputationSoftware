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

using namespace std;

/**Used 
   
 */
class CompositeBasisFunction
{
public:
  CompositeBasisFunction(vector<BasisFunction> _functions, ///Basis functions to use.
						 EigenInformation * myInformation, ///Eigeninformation. Ownership is not taken. May be deleted after constructor call.
						 const SpecificUnits * units ///Used when constructing.
						 ); ///Constructor. Ownership of the objects are NOT passed, and they are NOT used after the constructor.

  ~CompositeBasisFunction(); ///Destructor.

  ComplexDouble Eval(const double & x, ///x-value for evaluation.
					 uint pIndex ///Basis function index.
					 ); ///Evaluate the function using parameter x.

  ComplexDouble GetK(uint pIndex) const;

  ComplexDouble GetE(uint pIndex) const;

  uint GetNumberOfParameters() const; /// = number of eigenvalues, = number of basis states.

private:
  vector<BasisFunction> functions; ///Basis functions.
  vector<ComplexDouble> KValues; ///K-values.
  vector<ComplexDouble> Energies; ///Energies.
  vector<vector<ComplexDouble> > coefficients; ///Coefficients.
};

#endif
