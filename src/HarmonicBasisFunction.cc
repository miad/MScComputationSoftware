#include "HarmonicBasisFunction.hh"

HarmonicBasisFunction::HarmonicBasisFunction(vector<double> _xmin, vector<double> _omega, vector<Potential *> * _potentials, SpecificUnits * _units, uint precision)
  : xmin(_xmin), omega(_omega), potentials(_potentials), units(_units)
{
  if(units == NULL)
	throw RLException("Invalid units.");

  if(omega.size() != 2 || omega.at(0) <= 0 || omega.at(1) <= 0)
	{
	  throw RLException("Omega must be strictly greater than zero.");
	}
  if(potentials == NULL)
	{
	  throw RLException("Potential vector pointer cannot be null.");
	}
  for(vector<Potential*>::const_iterator it = potentials->begin(); it!=potentials->end(); ++it)
	{
	  if(*it == NULL)
		{
		  throw RLException("Potential vector cannot contain NULL elements.");
		}
	}
  if(xmin.size() != omega.size() || xmin.size() != potentials->size() )
	{
	  throw RLException("Invalid constructor vectors for HarmonicBasisFunction: code %d%d%d",xmin.size(),omega.size(),potentials->size());
	}


  for(uint i = 0; i<potentials->size(); ++i)
	GHpoints.push_back(HermiteRule::GetRule(precision, xmin.at(i), MASS*omega.at(i)/HBAR));

  InitFactorials();
}

void HarmonicBasisFunction::InitFactorials()
{
  factorials.push_back(1);
  for(uint i = 1; i<MAX_FACTORIALS; ++i)
	{
	  factorials.push_back(factorials[i-1]*i);
	}
  if(!DBL_EQUAL_MPREC(factorials[140], 1.34620124757175246058760738589E241) )
	{
	  throw RLException("Fatal error: factorial function doesn't perform as it should.");
	}
}

double HarmonicBasisFunction::GetEigenEnergy(uint n, uint pIndex) const
{
  return HBAR * omega.at(pIndex) * (0.5 + n) + potentials->at(pIndex)->Evaluate(xmin.at(pIndex));
}

HarmonicBasisFunction::~HarmonicBasisFunction()
{
  potentials = NULL;
  units = NULL;
}

long double HarmonicBasisFunction::Integrate(uint n1, uint n2, uint pIndex) const
{
  if(n1 > MAX_FACTORIALS || n2 > MAX_FACTORIALS)
	throw RLException("n > %d not implemented", MAX_FACTORIALS);


  long double value = 0;
  for(vector<pair<double, double> >::const_iterator it = GHpoints.at(pIndex).begin(); it!=GHpoints.at(pIndex).end(); ++it)
	{
	  long double WF1 = EvalNonExponentPart(n1, it->first, pIndex);
	  long double WF2 = EvalNonExponentPart(n2, it->first, pIndex); 
	  value += it->second * WF1 * WF2 * potentials->at(pIndex)->Evaluate(it->first);
	}

  return value;
}

 long double HarmonicBasisFunction::KineticTerm(uint n1, uint n2, uint pIndex) const
{
  if(n1 == n2)
	{
	  return HBAR * omega.at(pIndex) * 0.25 * (2.0 * n1 + 1.0);
	}
  else if(n1 == n2 + 2 || n1 + 2 == n2)
	{
	  uint n = MIN(n1, n2);
	  return -1* HBAR * omega.at(pIndex) * 0.25 * sqrt((n+1)*(n+2));
	}
  return 0;
}



long double HarmonicBasisFunction::DiffIntegrate(uint n1, uint n2, uint pIndex) const
{
  if(n1 > MAX_FACTORIALS || n2 > MAX_FACTORIALS)
	throw RLException("n > %d not implemented", MAX_FACTORIALS);


  long double value = 0;
  for(vector<pair<double, double> >::const_iterator it = GHpoints.at(pIndex).begin(); it!=GHpoints.at(pIndex).end(); ++it)
	{
	  long double WF1 = EvalNonExponentPart(n1, it->first, pIndex);
	  long double WF2 = EvalNonExponentPart(n2, it->first, pIndex); 
	  value += it->second * WF1 * WF2 *
		( potentials->at(pIndex)->Evaluate(it->first) - 
		  (  
		   0.5*MASS*pow(omega.at(pIndex),2.0) * pow((it->first - xmin.at(pIndex)),2.0) + potentials->at(pIndex)->Evaluate(xmin.at(pIndex)) 
		  ) 
		) ;
	}

  return value;
}

 long double HarmonicBasisFunction::NormalizationConstant(uint n, uint pIndex) const
{
  long double partA = 1.0/sqrt(pow(2.0L, n) * factorials.at(n));
  long double partB = pow(MASS * omega.at(pIndex) / (PI * HBAR), 0.25);
  return partA * partB;
}

long double HarmonicBasisFunction::EvalNonExponentPart(uint n, double x, uint pIndex) const
{
  long double part = HermiteEvaluator::HermiteH(n, sqrt(MASS*omega.at(pIndex)/HBAR) * (x-xmin.at(pIndex)));

  return part * NormalizationConstant(n, pIndex);
}


long double HarmonicBasisFunction::Eval(uint n, double x, uint pIndex) const
{

  long double partD = exp(-MASS*omega.at(pIndex)*pow(x-xmin.at(pIndex),2.0)/(2.0*HBAR));
  return EvalNonExponentPart(n, x, pIndex) * partD;
}

double HarmonicBasisFunction::GetXmin(uint pIndex ) const
{
  return xmin.at(pIndex);
}

double HarmonicBasisFunction::GetOmega(uint pIndex) const
{
  return omega.at(pIndex);
}
