#include "ParametrizedPotential.hh"

ParametrizedPotential::ParametrizedPotential(string function, vector<pair<string, double> > parameters, double minX, double maxX)
  :minX(minX), maxX(maxX)
{
  SetPrecision(100);
  for(vector<pair<string, double> >::const_iterator it = parameters.begin(); it!=parameters.end(); ++it)
	{
	  fp.AddConstant(it->first, it->second);
	}
  fp.Parse(function, "x");

  fp.Optimize();
}

ParametrizedPotential::ParametrizedPotential(const ParametrizedPotential & other)
  :fp(other.fp), legendreRule(other.legendreRule), minX(other.minX), maxX(other.maxX), precision(other.precision)
{
  fp.ForceDeepCopy();
}

ParametrizedPotential::~ParametrizedPotential()
{

}

double ParametrizedPotential::GetMinX() const
{
  return minX;
}

double ParametrizedPotential::GetMaxX() const
{
  return maxX;
}

double ParametrizedPotential::Evaluate(double x) const
{
  return const_cast<FunctionParser*>(&fp)->Eval(&x);
}


ComplexDouble ParametrizedPotential::BasisIntegrate(BasisFunction & b1, BasisFunction & b2, ComplexDouble & k1, ComplexDouble & k2)
{
  ComplexDouble value = ComplexDouble(0.0, 0.0);

  for(vector<pair<double, double> >::const_iterator it = legendreRule.begin(); it!=legendreRule.end(); ++it)
	{
	  value += it->second * fp.Eval(&it->first) * b1.Eval(it->first, k1) * b2.Eval(it->first, k2);
	}

  return value;
}


void ParametrizedPotential::SetPrecision(unsigned long value)
{
  precision = value;
  legendreRule = LegendreRule::GetRule(precision, minX, maxX);
}


unsigned long ParametrizedPotential::GetPrecision() const
{
  return precision;
}

vector<pair<double, double> > ParametrizedPotential::GetPlottingPoints() const
{
  vector<pair<double, double> > toReturn;
  
  toReturn.push_back(make_pair(minX, 0.0));
  for(vector<pair<double, double> >::const_iterator it = legendreRule.begin(); it!=legendreRule.end(); ++it)
	{
	  toReturn.push_back(make_pair(it->first, const_cast<FunctionParser*>(&fp)->Eval(&it->first)));
	}
  toReturn.push_back(make_pair(maxX, 0.0));

  return toReturn;
}


vector<pair<double, double> > ParametrizedPotential::GetPrecisionPoints() const
{
  vector<pair<double, double> > toReturn;
  
  for(vector<pair<double, double> >::const_iterator it = legendreRule.begin(); it!=legendreRule.end(); ++it)
	{
	  toReturn.push_back(make_pair(it->first, const_cast<FunctionParser*>(&fp)->Eval(&it->first)));
	}

  return toReturn;
}


PotentialType ParametrizedPotential::GetType() const
{
  return TypeParametrized;
}