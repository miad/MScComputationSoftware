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

using namespace std;



class BasisFunction
{
public:
  BasisFunction(string _name ///Name of the basis function. May be arbitrary function.
				);///Constructor.
  ComplexDouble Eval(const ComplexDouble & x, 
					 const ComplexDouble & k = 1.0
					 );
  const char * GetName() const;
  
  void ForceDeepCopy(); ///Forces a deep copy for the underlying FunctionParser object.

  //ComplexDouble GetPreFactor() const;
  


private:
  FunctionParser_cd fp;
  //ComplexDouble factor; ///Indicates sign and possibly a factor I.
  //ComplexDouble preFactor; ///A factor in front of the basis element. 
  short type; ///Indicates the type (sin, cos, exp, ...)
  string name;
};

#endif
