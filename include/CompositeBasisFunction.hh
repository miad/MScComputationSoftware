#ifndef CompositeBasisFunction_hh
#define CompositeBasisFunction_hh 1

#include "Globals.hpp"
#include "RLException.hh"
#include "fparser.hh"
#include <stdio.h>
#include <ctype.h>
#include <cstring>
#include <map>
#include <iostream>
#include "HermiteEvaluator.hh"
#include "BasisFunction.hh"

using namespace std;

/**Used 
   
 */
class CompositeBasisFunction
{
public:
  CompositeBasisFunction(vector<BasisFunction> * _functions, ///Basis functions to use.
						 vector<ComplexDouble> * _parameters, ///Essentially the eigenvalues (converted from energy to k-values of course).
						 vector<vector<ComplexDouble> > * _coefficients ///The coefficients. Basically 1-particle eigenvectors.
						 ); ///Constructor. Ownership of the objects are NOT passed, but they may of course NOT be deleted before the deletion of this object.

  ~CompositeBasisFunction(); ///Destructor.

  ComplexDouble Eval(const double & x, ///x-value for evaluation.
					 uint pIndex ///Basis function index.
					 ); ///Evaluate the function using parameter x.

  uint GetNumberOfParameters() const; /// = number of eigenvalues, = number of basis states.
  
private:
  vector<BasisFunction> * functions; ///Functions.
  vector<ComplexDouble> * parameters; ///Parameters.
  vector<vector<ComplexDouble> > * coefficients; ///Coefficients.
};

#endif
