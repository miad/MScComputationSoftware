#include "Potential.hh"

Potential::Potential()
  :minX(0.), maxX(0.)
{

}

Potential::~Potential()
{

}

ComplexDouble Potential::FastExpIntegrate(ComplexDouble exponentVal)
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

void Potential::AddValue(Interval toAdd)
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

void Potential::AddValue(double x1, double x2, double y)
{
  AddValue(Interval(x1, x2, y));
}

double Potential::Evaluate(double x) const
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

list<Interval> Potential::GetPotentialPoints() const
{
  return PotentialPoints;
}

double Potential::GetMinX() const
{
  return minX;
}

double Potential::GetMaxX() const
{
  return maxX;
}

unsigned int Potential::GetNumberOfValues() const
{
  return PotentialPoints.size();
}
