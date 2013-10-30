#ifndef Interval_hh
#define Interval_hh 1
#include "RLException.hh"
#include "Globals.hpp"
#include <algorithm>

using namespace std;

class Interval
{
public:
  Interval(double _x1, double _x2, double _y);
  ~Interval();
  bool Overlaps(const Interval & toCheck) const;
  double x1, x2, y;
};

#endif
