#include "PiecewiseConstantPotential.hh"

PiecewiseConstantPotential::PiecewiseConstantPotential()
  :minX(0.), maxX(0.), precision(100)
{

}

PiecewiseConstantPotential::PiecewiseConstantPotential(string filename)
{
  FILE * inputFile = fopen(filename.c_str(), "r");
  char * line = NULL; ///stores line data when read from file.
  size_t len = 0; ///is ignored. 
  ssize_t numberOfReadCharacters; ///length of read text

  if(inputFile == NULL)
	throw RLException("Could not open potential input file '%s'.", filename.c_str());

  while ((numberOfReadCharacters = getline(&line, &len, inputFile)) != -1) 
	{
	  ///check if we should ignore the line.
	  int ptr = 0;
	  while(ptr < numberOfReadCharacters && line[ptr] == ' ')
		++ptr;
	  if(ptr == numberOfReadCharacters || line[ptr] == '#' || line[ptr] == '\n')
		continue;

	  double x1, x2, y;
	  if(sscanf(line, "%lf %lf %lf",&x1, &x2, &y) == EOF)
		{
		  throw RLException("Suspicious line in potential file '%s': '%s'", filename.c_str(), line);
		}
	  AddValue(x1, x2, y);
	}

  if (line)
	free(line);
  fclose(inputFile);
  if(PotentialPoints.empty())
	throw RLException("The specified potential file '%s' does not contain any potential data.", filename.c_str());
}

PiecewiseConstantPotential::~PiecewiseConstantPotential()
{

}

ComplexDouble PiecewiseConstantPotential::FastExpIntegrate(ComplexDouble exponentVal) const
{
  ComplexDouble value = 0;
  if(abs(exponentVal) < EPS) ///If the exponent is zero, the integrand is computed differently.
	{
	  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
		{
		  value += it->y*(it->x2 - it->x1);
		}
	  return value;
	}

  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{
	  value+=it->y*(exp(ComplexDouble(0, 1)*exponentVal*it->x2) - 
					exp(ComplexDouble(0, 1)*exponentVal*it->x1) );
	}
  value /= ComplexDouble(0, 1)*exponentVal;

  return value;
}

ComplexDouble PiecewiseConstantPotential::FastCosIntegrate(ComplexDouble k1, ComplexDouble k2) const
{  
  ComplexDouble value = 0;
  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{
	  if( DBL_EQUAL(k1, k2) && DBL_EQUAL(k1, 0.)) ///If the exponent is zero, the integrand is computed differently.
		{
		  value += it->y*(it->x2 - it->x1);
		}
	  else if( DBL_EQUAL(k1, k2) || DBL_EQUAL(k1, -k2) )
		{
		  value += it->y * ( (it->x2 - it->x1)*0.5 +
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else
		{
		  value += it->y * 0.5*(
								(sin((k1-k2)*it->x2)-sin((k1-k2)*it->x1))/(k1-k2) +
								(sin((k1+k2)*it->x2)-sin((k1+k2)*it->x1))/(k1+k2)
								);
		}
	}
  return value;
}

ComplexDouble PiecewiseConstantPotential::FastSinIntegrate(ComplexDouble k1, ComplexDouble k2) const
{

  ComplexDouble value = 0;
  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{	  
	  if( DBL_EQUAL(k1, k2) && DBL_EQUAL(k1, 0.)) ///If the exponent is zero, the integrand is computed differently.
		{
		  value = 0.;
		  return value;
		}
	  else if( DBL_EQUAL(k1, k2))
		{

		  value += it->y * ( (it->x2 - it->x1)*0.5 -
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else if( DBL_EQUAL(k1, -k2))
		{

		  value += -1.*it->y * ( (it->x2 - it->x1)*0.5 -
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else
		{

		  value += it->y * 0.5*(
								(sin((k1-k2)*it->x2)-sin((k1-k2)*it->x1))/(k1-k2) -
								(sin((k1+k2)*it->x2)-sin((k1+k2)*it->x1))/(k1+k2)
								);
		}
	}
  return value;
}


void PiecewiseConstantPotential::AddValue(Interval toAdd)
{
  foric(list<Interval>, PotentialPoints, it)
	{
	  if(it->Overlaps(toAdd))
		{
		  throw RLException("Attempted to add an overlapping interval.");
		}
	}
  minX = PotentialPoints.empty()?toAdd.x1:MIN(minX, toAdd.x1);
  maxX = PotentialPoints.empty()?toAdd.x2:MAX(maxX, toAdd.x2);

  PotentialPoints.push_back(toAdd);
}

void PiecewiseConstantPotential::AddValue(double x1, double x2, double y)
{
  AddValue(Interval(x1, x2, y));
}

double PiecewiseConstantPotential::Evaluate(double x) const
{
  foric(list<Interval>, PotentialPoints, it)
	{
	  if( it->x1 <= x && it->x2 > x )
		{
		  return it->y;
		}
	}
  return 0;
}

void PiecewiseConstantPotential::RecomputeLegendreRules()
{
  legendreRules.clear();
  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!= PotentialPoints.end(); ++it)
	{
	  unsigned int pointsOnCurve = 1 + (int)(precision*((it->x2-it->x1)/(maxX-minX)));
	  legendreRules.push_back(LegendreRule::GetRule(pointsOnCurve, it->x1, it->x2));
	}
}

list<Interval> PiecewiseConstantPotential::GetPotentialPoints() const
{
  return PotentialPoints;
}

double PiecewiseConstantPotential::GetMinX() const
{
  return minX;
}

double PiecewiseConstantPotential::GetMaxX() const
{
  return maxX;
}

unsigned int PiecewiseConstantPotential::GetNumberOfValues() const
{
  return PotentialPoints.size();
}


void PiecewiseConstantPotential::Clear() 
{
  PotentialPoints.clear();
  minX = 1E99;
  maxX = -1E99;
}


ComplexDouble PiecewiseConstantPotential::BasisIntegrate(BasisFunction & b1, BasisFunction & b2, ComplexDouble & k1, ComplexDouble & k2)
{
  if(legendreRules.size() != PotentialPoints.size() )
	{
	  throw RLException("The legendre rules for the piecewise potential were not computed. Call RecomputeLegendreRules before trying to evaluate.");
	}

  ///Plain ol' Gauss-Legendre integration.
  list<Interval>::const_iterator iPot = PotentialPoints.begin();
  list<vector<pair<double, double> > >::const_iterator iLeg = legendreRules.begin();
  ComplexDouble value = ComplexDouble(0.0,0.0);
  while(iPot != PotentialPoints.end())
  {
	for(vector<pair<double, double> >::const_iterator it = iLeg->begin(); it != iLeg->end(); ++it)
	  {
		value += it->second * iPot->y * b1.Eval(it->first, k1) * b2.Eval(it->first, k2);
	  }
	++iPot; ++iLeg;
  }
  return value;
}


void PiecewiseConstantPotential::SetPrecision(unsigned long value)
{
  precision = value;
}


unsigned long PiecewiseConstantPotential::GetPrecision() const
{
  return precision;
}

vector<pair<double, double> > PiecewiseConstantPotential::GetPlottingPoints() const
{
  vector<pair<double, double> > toReturn;

  list<Interval> localIntervals = PotentialPoints; ///We may not modify inside a const function, so copy.
  localIntervals.sort();
  
  toReturn.push_back(make_pair(localIntervals.front().x1, 0.0));

  for(list<Interval>::const_iterator it = localIntervals.begin(); it!=localIntervals.end(); ++it)
	{
	  ///The small offset EPS is to make the function injective.
	  toReturn.push_back(make_pair(it->x1 + EPS, it->y));
	  toReturn.push_back(make_pair(it->x2 - EPS,it->y));
	}
  toReturn.push_back(make_pair(localIntervals.back().x2, 0.0));

  return toReturn;
}


vector<pair<double, double> > PiecewiseConstantPotential::GetPrecisionPoints() const
{
  vector<pair<double, double> > toReturn;
  list<Interval>::const_iterator iPtr = PotentialPoints.begin();
  for(list<vector<pair<double, double> > >::const_iterator it = legendreRules.begin(); it!=legendreRules.end(); ++it)
	{
	  for(vector<pair<double, double> >::const_iterator ip = it->begin(); ip != it->end(); ++ip)
		{
		  toReturn.push_back(make_pair(ip->first, iPtr->y));
		}
	  ++iPtr;
	}
  return toReturn;

  return toReturn;
}



PotentialType PiecewiseConstantPotential::GetType() const
{
  return TypePiecewiseConstant;
}
