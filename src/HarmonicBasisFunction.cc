#include "HarmonicBasisFunction.hh"

HarmonicBasisFunction::HarmonicBasisFunction(double _xmin, double _omega, Potential * _pot, SpecificUnits * _units)
  : xmin(_xmin), omega(_omega), pot(_pot), units(_units)
{
  if(pot == NULL)
	throw RLException("Invalid potential.");
  if(units == NULL)
	throw RLException("Invalid units.");
  if(omega <= 0)
	{
	  throw RLException("Omega must be strictly greater than zero.");
	}
  
  GHpoints = HermiteRule::GetRule(pot->GetPrecision(), xmin, MASS*omega/HBAR);

  //GHpoints = LegendreRule::GetRule(1000, -20, 20); ///TEMP

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

double HarmonicBasisFunction::GetEigenEnergy(uint n) const
{
  return HBAR * omega * (0.5 + n) + pot->Evaluate(xmin);
}

HarmonicBasisFunction::~HarmonicBasisFunction()
{
  pot = NULL;
  units = NULL;
}

long double HarmonicBasisFunction::Integrate(uint n1, uint n2) const
{
  if(n1 > MAX_FACTORIALS || n2 > MAX_FACTORIALS)
	throw RLException("n > %d not implemented", MAX_FACTORIALS);


  long double value = 0;
  for(vector<pair<double, double> >::const_iterator it = GHpoints.begin(); it!=GHpoints.end(); ++it)
	{
	  long double WF1 = EvalNonExponentPart(n1, it->first);
	  long double WF2 = EvalNonExponentPart(n2, it->first); 
	  value += it->second * WF1 * WF2 * pot->Evaluate(it->first);
	}

  return value;
}

long double HarmonicBasisFunction::KineticTerm(uint n1, uint n2) const
{
  if(n1 == n2)
	{
	  return HBAR * omega * 0.25 * (2.0 * n1 + 1.0);
	}
  else if(n1 == n2 + 2 || n1 + 2 == n2)
	{
	  uint n = MIN(n1, n2);
	  return -1* HBAR * omega * 0.25 * sqrt((n+1)*(n+2));
	}
  return 0;
}



long double HarmonicBasisFunction::DiffIntegrate(uint n1, uint n2) const
{
  if(n1 > MAX_FACTORIALS || n2 > MAX_FACTORIALS)
	throw RLException("n > %d not implemented", MAX_FACTORIALS);


  long double value = 0;
  for(vector<pair<double, double> >::const_iterator it = GHpoints.begin(); it!=GHpoints.end(); ++it)
	{
	  long double WF1 = EvalNonExponentPart(n1, it->first);
	  long double WF2 = EvalNonExponentPart(n2, it->first); 
	  value += it->second * WF1 * WF2 *
		( pot->Evaluate(it->first) - 
		  (  
		   0.5*MASS*pow(omega,2.0) * pow((it->first - xmin),2.0) + pot->Evaluate(xmin) 
		  ) 
		) ;
	}

  return value;
}

long double HarmonicBasisFunction::NormalizationConstant(uint n) const
{
  long double partA = 1.0/sqrt(pow(2.0L, n) * factorials.at(n));
  long double partB = pow(MASS * omega / (PI * HBAR), 0.25);
  return partA * partB;
}

long double HarmonicBasisFunction::EvalNonExponentPart(uint n, double x) const
{
  long double part = HermiteEvaluator::HermiteH(n, sqrt(MASS*omega/HBAR) * (x-xmin));

  return part * NormalizationConstant(n);
}


long double HarmonicBasisFunction::Eval(uint n, double x) const
{

  long double partD = exp(-MASS*omega*pow(x-xmin,2.0)/(2.0*HBAR));
  return EvalNonExponentPart(n, x) * partD;
}

double HarmonicBasisFunction::GetXmin() const
{
  return xmin;
}

double HarmonicBasisFunction::GetOmega() const
{
  return omega;
}