#ifndef ParametrizedPotential_hh
#define ParametrizedPotential_hh 1

#include <list>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "Interval.hh"
#include "RLException.hh"
#include "Globals.hpp"
#include "Potential.hh"
#include "LegendreRule.hh"

using namespace std;

/** 
	Stepwise potential.
	Todo: implement arbitrary function between the points in the potential.
	Will be needed later on.
 */
class ParametrizedPotential
  : public Potential
{
public:
  ParametrizedPotential(string function, /// The function to use.
						vector<pair<string, double> > parameters, /// The parameters to use. 
						double minX, /// Minimum x value.
						double maxX /// Maximum x value.
						); ///Constructor. \todo Implement reading potential from file.

  ParametrizedPotential(string function, /// The function to use.
						vector<pair<string, double> > parameters, /// The parameters to use. 
						const char * minXFunction, ///Min x function.
						const char * maxXFunction ///Max x function.
						); ///Constructor. \todo Implement reading potential from file.

  ParametrizedPotential(const ParametrizedPotential &other ///Other potential
						); ///Copy constructor.

  ~ParametrizedPotential(); ///Destructor.

  double Evaluate(double x ///The point to evaluate the potential value at.
				  ) const; ///Evaluate a specific point in this potential.

  double GetMinX() const; ///Returns the minimum x value for which the potential is nonzero.
  double GetMaxX() const; ///Returns the maximum x value for which the potential is nonzero.

  unsigned int GetNumberOfValues() const;

  list<Interval> GetPotentialPoints() const;

  ComplexDouble BasisIntegrate(BasisFunction & b1, ///First basis function.
							   BasisFunction & b2, ///Second basis function.
							   ComplexDouble & k1, ///First k-value.
							   ComplexDouble & k2 ///Second k-value.
							   ) ; ///Integrate the basis functions over the potential to create V(k1, k2). 
  
  
  unsigned long GetPrecision() const; ///Returns the precision.

  void SetPrecision(unsigned long value ///The value of the new precision.
					); ///Set the precision to the value given by 'value'.

  vector<pair<double, double> > GetPlottingPoints() const; ///Returns suitable points for plotting.

  vector<pair<double, double> > GetPrecisionPoints() const; ///Returns points suitable for plotting together with the potential, indicating precision.

  PotentialType GetType() const;


private:
  FunctionParser fp; ///Function parser used to evaluate the potential in a specific point.

  vector<pair<double, double> > legendreRule;

  double minX, maxX; ///Minimum and maximum x-value.
  unsigned long precision; ///Precision. Here: number of GL points on each curve.
};



#endif
