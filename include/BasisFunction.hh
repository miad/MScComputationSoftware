#ifndef BasisFunction_hh
#define BasisFunction_hh 1

#include "Globals.hpp"
#include "RLException.hh"
#include "fparser.hh"
#include <stdio.h>
#include <ctype.h>
#include <cstring>
#include <map>
#include <iostream>
#include "HermiteEvaluator.hh"

using namespace std;



class BasisFunction
{
public:
  BasisFunction(string _name ///Name of the basis function. May be arbitrary function.
				);///Constructor.

  BasisFunction(const BasisFunction & other ///To copy from.
				); ///Copy constructor.


  ComplexDouble Eval(const ComplexDouble & x, ///The x-value for evaluation.
					 const ComplexDouble & k = 1.0 ///The k-value for evaluation.
					 ); ///Evaluate the function with the parameters x and k.

  const char * GetName() const; ///Returns a string representing the function.
  
private:
  FunctionParser_cd fp; ///Function parser used to parse.

  string name; ///Contains a string representing the type of basis function. This is the same as was sent to the constructor.
};

#endif
