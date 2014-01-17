#include "CompositeBasisFunction.hh"

CompositeBasisFunction::CompositeBasisFunction(vector<BasisFunction> * _functions, vector<ComplexDouble> *_parameters, vector<vector<ComplexDouble> > * _coefficients)
  : functions(_functions), parameters(_parameters), coefficients(_coefficients)
{
  if(functions == NULL || parameters == NULL || coefficients == NULL)
	{
	  throw RLException("Not allowed to use NULL argument for a CompositeBasisFunction constructor.");
	}
  if(coefficients->size() != parameters->size() || coefficients->size() % functions->size() != 0)
	{
	  throw RLException("Invalid input array size for composite basis function.");
	}
}

CompositeBasisFunction::~CompositeBasisFunction()
{
  ///Don't delete since we don't have ownership.
  functions = NULL;
  parameters = NULL;
  coefficients = NULL;
}


ComplexDouble CompositeBasisFunction::Eval(const double & x, uint pIndex)
{
  if(pIndex >= parameters->size())
	{
	  throw RLException("Invalid p-index.");
	}

  uint bConvert = coefficients->size() / functions->size();
  
  ComplexDouble sum = 0.0;

  for(uint i = 0; i<coefficients->at(pIndex).size(); ++i)
	{
	  uint basisPointer = i/bConvert;
	  sum += coefficients->at(pIndex).at(i) * functions->at(basisPointer).Eval(x, parameters->at(pIndex));
	}

  return sum;
}

uint CompositeBasisFunction::GetNumberOfParameters() const
{
  return parameters->size();
}
