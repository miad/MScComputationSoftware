#include "Interval.hh"
Interval::Interval(double _x1, double _x2, double _y)
  :x1(_x1), x2(_x2), y(_y)
{
  if( x1 > x2 )
	swap(x1, x2);
  if( x1 == x2 )
	throw basisException("Invalid interval: length was zero");
}

Interval::~Interval()
{
  
}

bool Interval::Overlaps(const Interval & toCheck) const
{
  return ( toCheck.x1 > x1 && toCheck.x1 < x2 ) || ( toCheck.x2 > x1 && toCheck.x2 < x2 );
}
