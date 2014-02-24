#include "CompositeBasisFunction.hh"

CompositeBasisFunction::CompositeBasisFunction(vector<BasisFunction> _functions, EigenInformation * _myInformation, const ParametrizedCurve * _KCurve)
  : functions(_functions), myHarmonicBasisFunction(NULL), myInformation(_myInformation), KCurve(_KCurve)
{
  if(myInformation == NULL || KCurve == NULL || functions.empty())
	{
	  throw RLException("CompositeBasisFunction: tried to initialize with too empty input parameters.");
	}

  if(myInformation->Eigenvalues.size() % functions.size() != 0)
	{
	  throw RLException("CompositeBasisFunction: Invalid input dimensions.");
	}
}

CompositeBasisFunction::CompositeBasisFunction(const HarmonicBasisFunction * _myHarmonicBasisFunction, EigenInformation * _myInformation, uint _curveIndex)
  : myHarmonicBasisFunction(_myHarmonicBasisFunction), myInformation(_myInformation), KCurve(NULL), curveIndex(_curveIndex)
{
  if(myInformation == NULL || myHarmonicBasisFunction == NULL)
	{
	  throw RLException("CompositeBasisFunction: tried to call with NULL EigenInformation or HarmonicBasisFunction object.");
	}
}


CompositeBasisFunction::~CompositeBasisFunction()
{
  ///Don't delete since we don't have ownership.
}

ComplexDouble CompositeBasisFunction::GetE(uint pIndex) const
{
  return myInformation->Eigenvalues.at(pIndex);
}


ComplexDouble CompositeBasisFunction::Eval(const double & x, uint pIndex)
{
  if(pIndex >= myInformation->Eigenvectors.size())
	{
	  throw RLException("Invalid p-index: %d but %d is max", pIndex, myInformation->Eigenvectors.size());
	}

  if(myHarmonicBasisFunction != NULL)
	return HarmonicEval(x, pIndex);
  return KEval(x, pIndex);
  
}

ComplexDouble CompositeBasisFunction::KEval(const double & x, uint pIndex)
{
  uint N = KCurve->GetTotalNumberOfGLPoints();
  ComplexDouble sum = 0;

  for(uint i = 0; i<myInformation->Eigenvectors.at(pIndex).size(); ++i)
	{
	  ///By convention (when creating the Hamiltonian)
	  uint curvePointer = i % N;
	  uint basisPointer = i / N;
	  ComplexDouble kVal = KCurve->GetRuleValue(curvePointer);
	  ComplexDouble kWeight = KCurve->GetRuleWeight(curvePointer);
	  
	  sum += myInformation->Eigenvectors.at(pIndex).at(i) * 
		sqrt(kWeight) * functions.at(basisPointer).Eval(x, kVal);
	}
  return sum;
}



double CompositeBasisFunction::HarmonicEval(const double & x, uint pIndex) const
{
  long double sum = 0;
  for(uint n = 0; n<myInformation->Eigenvectors.at(pIndex).size(); ++n)
	{
	  ComplexDouble coefficient = myInformation->Eigenvectors.at(pIndex).at(n);
	  if(imag(coefficient) > 1E-9)
		throw RLException("Complex coefficient (imag part %+13.10f ) for harmonic oscillator basis: something is wrong here.", imag(coefficient));

	  sum += real(coefficient) * myHarmonicBasisFunction->Eval(n, x, curveIndex);
	}
  return (double)sum;
}

uint CompositeBasisFunction::GetSize() const
{
  return myInformation->Eigenvalues.size();
}


bool CompositeBasisFunction::IsHarmonic() const
{
  return myHarmonicBasisFunction != NULL;
}

pair<uint, ComplexDouble> CompositeBasisFunction::GetDominatingVectorPart(uint eigenVector) const
{
  pair<uint, ComplexDouble> best = make_pair(0, 0.0);
  for(uint i = 0; i<myInformation->Eigenvectors.at(eigenVector).size(); ++i)
	{
	  if(abs(myInformation->Eigenvectors.at(eigenVector).at(i)) > abs(best.second))
		{
		  best.second = myInformation->Eigenvectors.at(eigenVector).at(i);
		  best.first = i;
		}
	}
  return best;
}
