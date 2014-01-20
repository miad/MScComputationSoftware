#include "CompositeBasisFunction.hh"

CompositeBasisFunction::CompositeBasisFunction(vector<BasisFunction> _functions, EigenInformation * myInformation, const SpecificUnits * units)
  : functions(_functions)
{
  if(myInformation == NULL || units == NULL)
	{
	  throw RLException("CompositeBasisFunction: tried to call with NULL EigenInformation or units object.");
	}

  Energies = vector<ComplexDouble>(myInformation->Eigenvalues);
  for(vector<ComplexDouble>::const_iterator it = Energies.begin(); it!=Energies.end(); ++it)
	{
	  KValues.push_back(units->EnergyToKValue(*it));
	}

  coefficients = vector<vector<ComplexDouble> >(myInformation->Eigenvectors);


if(coefficients.size() != KValues.size() || coefficients.size() % functions.size() != 0)
  {
	throw RLException("Invalid input array size for composite basis function.");
  }
}

CompositeBasisFunction::~CompositeBasisFunction()
{
  ///Don't delete since we don't have ownership.
}

ComplexDouble CompositeBasisFunction::GetE(uint pIndex) const
{
  return Energies.at(pIndex);
}

ComplexDouble CompositeBasisFunction::GetK(uint pIndex) const
{
  return KValues.at(pIndex);
}


ComplexDouble CompositeBasisFunction::Eval(const double & x, uint pIndex)
{
  if(pIndex >= KValues.size())
	{
	  throw RLException("Invalid p-index.");
	}
  
  uint bConvert = coefficients.size() / functions.size();
  
  ComplexDouble sum = 0.0;
  
  for(uint i = 0; i<coefficients.at(pIndex).size(); ++i)
	{
	  uint basisPointer = i/bConvert;
	  sum += coefficients.at(pIndex).at(i) * functions.at(basisPointer).Eval(x, KValues.at(pIndex));
	}


  return sum;
}

uint CompositeBasisFunction::GetNumberOfParameters() const
{
  return KValues.size();
}
