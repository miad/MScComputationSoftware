#ifndef BasisFunction_hh
#define BasisFunction_hh 1

#include "Globals.hpp"
#include "RLException.hh"
#include <stdio.h>
#include <ctype.h>
#include <cstring>

using namespace std;



class BasisFunction
{
public:
  BasisFunction(const char * name ///Name of the basis function. May be one of the following:
				/*    "sin"
				 *    "cos"
				      "exp"
					  
					  also, the fourth and fifth characters may be a plus sign (implicit by omission), a minus sign and/or an i.
					  Normal choices: 
					  expi+ expi-
					  sin cos
				 */
				);///Constructor.
  ComplexDouble Eval(ComplexDouble x);
  const char * GetName() const;

private:
  ComplexDouble factor; ///Indicates sign and possibly a factor I.
  short type; ///Indicates the type (sin, cos, exp, ...)
  char name[6]; ///Used when returning name.
};

#endif
