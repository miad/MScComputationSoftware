#include "Potential.hh"

Potential::Potential()
{
  minX = 1E99; maxX = -1E99;
  InitializeDefault();
}

Potential::~Potential()
{

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
  minX = MIN(minX, toAdd.x1);
  maxX = MAX(maxX, toAdd.x2);

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
	  if( it->x1 < x && it->x2 > x )
		{
		  return it->y;
		}
	}
  return 0;
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

void Potential::InitializeDefault()
{
  //Define the shape of the potential.
  //NOTE: Patches must be in order and non-overlapping.
  AddValue(-1.5, -0.5, 5);
  AddValue(-0.5, 0.5, -10);
  AddValue(0.5, 1.5, 5);
}

