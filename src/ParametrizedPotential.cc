#include "ParametrizedPotential.hh"

ParametrizedPotential::ParametrizedPotential(string function, vector<pair<string, double> > parameters, double minX, double maxX)
  :minX(minX), maxX(maxX)
{
  SetPrecision(100);
  for(vector<pair<string, double> >::const_iterator it = parameters.begin(); it!=parameters.end(); ++it)
	{
	  if(!fp.AddConstant(it->first, it->second))
		{
		  throw RLException("Could not add parameter '%s' with value '%d' to ParametrizedPotential function.", it->first.c_str(), it->second);
		}
	}
  int retVal = fp.Parse(function, "x");
  if(retVal != -1)
	{
	  throw RLException("Error in parsing the function '%s' at character %d.", function.c_str(), retVal);
	}

  fp.Optimize();
}

ParametrizedPotential::ParametrizedPotential(string function, vector<pair<string, double> > parameters, const char * minXFunction, const char * maxXFunction)
{
  SetPrecision(100);
  for(vector<pair<string, double> >::const_iterator it = parameters.begin(); it!=parameters.end(); ++it)
	{
	  fp.AddConstant(it->first, it->second);
	}

  int retVal1 = fp.Parse(minXFunction, "");
  if(retVal1 != -1)
	{
	  char buffer[10000];
	  memset(buffer, ' ', sizeof(char)*9995);
	  buffer[retVal1 +1] = '\0';
	  throw RLException("Error in parsing the function \n'%s' at character %d\n%s%c\n%s", 
						function.c_str(), 
						retVal1,
						buffer,
						'#',
						fp.ErrorMsg());
	}

  minX = fp.Eval(NULL);
  int retVal2 = fp.Parse(maxXFunction, "");
  if(retVal2 != -1)
	{
	  char buffer[10000];
	  memset(buffer, ' ', sizeof(char)*9995);
	  buffer[retVal2 +1] = '\0';
	  throw RLException("Error in parsing the function \n'%s' at character %d\n%s%c\n%s", 
						function.c_str(), 
						retVal2,
						buffer,
						'#',
						fp.ErrorMsg());

	}

  maxX = fp.Eval(NULL);


  int retVal3 = fp.Parse(function, "x");
  if(retVal3 != -1)
	{
	  char buffer[10000];
	  memset(buffer, ' ', sizeof(char)*9995);
	  buffer[retVal3 +1] = '\0';
	  throw RLException("Error in parsing the function \n'%s' at character %d\n%s%c\n%s", 
						function.c_str(), 
						retVal3,
						buffer,
						'#',
						fp.ErrorMsg());

	}
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
  if(x < minX || x > maxX)
	return 0;
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
