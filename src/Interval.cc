#include "Interval.hh"
Interval::Interval(double _x1, double _x2, double _y)
  :x1(_x1), x2(_x2), y(_y)
{
  if( x1 > x2 )
	swap(x1, x2);
  if( x1 == x2 )
	throw RLException("Invalid interval: length was zero");
}

Interval::~Interval()
{
  
}

bool Interval::Overlaps(const Interval & toCheck) const
{
  return ( toCheck.x1 > x1 && toCheck.x1 < x2 ) || ( toCheck.x2 > x1 && toCheck.x2 < x2 );
}

const bool Interval::operator<(const Interval &other) const
{
  if( !DBL_EQUAL(x1, other.x1))
	return x1 < other.x1;
  if( !DBL_EQUAL(x2, other.x2))
	return x2 < other.x2;
  if( !DBL_EQUAL(y, other.y))
	return y < other.y;
  return false;
}

const bool Interval::operator==(const Interval &other) const
{
  return ! ((other < *this) || (*this < other) );
}

const bool Interval::operator>(const Interval &other) const
{
  return other < *this;
}

const bool Interval::operator>=(const Interval &other) const
{
  return ((*this > other ) || (*this == other));
}

const bool Interval::operator<=(const Interval &other) const
{
  return ((*this < other) || (*this > other));
}

const bool Interval::operator!=(const Interval &other) const
{
  return !(*this==other);
}
