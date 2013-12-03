#include "BasisFunction.hh"

BasisFunction::BasisFunction(string _name)
  :type(-1), name(_name)
{

  fp.Parse(name, "x,k");
  fp.Optimize();
}

BasisFunction::BasisFunction(const BasisFunction & other)
{
  fp = other.fp;
  type = other.type;
  name = other.name;
  fp.ForceDeepCopy(); ///Essential for thread safety.
}

ComplexDouble BasisFunction::Eval(const ComplexDouble & x, const ComplexDouble & k)
{
  ComplexDouble toEval[2]; toEval[0] = x; toEval[1] = k;
  ComplexDouble toReturn = fp.Eval(&toEval[0]);
  int err = fp.EvalError();
  if(err)
	{
	  throw RLException("Function evaluation error: %d\n", err);
	}
  return toReturn;
}

const char * BasisFunction::GetName() const
{
  return name.c_str();
}
