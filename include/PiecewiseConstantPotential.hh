#ifndef PiecewiseConstantPotential_hh
#define PiecewiseConstantPotential_hh 1

#include <list>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include "Interval.hh"
#include "RLException.hh"
#include "Globals.hpp"
#include "Potential.hh"

using namespace std;

/** 
	Stepwise potential.
	Todo: implement arbitrary function between the points in the potential.
	Will be needed later on.
 */
class PiecewiseConstantPotential
  : public Potential
{
public:
  PiecewiseConstantPotential(); ///Constructor. \todo Implement reading potential from file.
  PiecewiseConstantPotential(string filename ///name of file.
			); ///Construct a potential from a space-separated file. The file must consist of lines with 3 rows and one double precision number in each of these. Lines starting with # are ignored.
  ~PiecewiseConstantPotential(); ///Destructor.
  void AddValue(double x1, ///The x1 value to add.
				double x2, ///The x2 value to add.
				double y ///The y value to add.
				); ///Adds a value to the stepwise potential.

  ComplexDouble FastExpIntegrate(ComplexDouble exponentVal ///The exponent in the fast integration.
								 ) const; ///Fast integration of e^{i*exponentVal*x} over the potential. Returns the integral over the entire potential.

  ComplexDouble FastCosIntegrate(ComplexDouble k1, ///Coefficient 1
								 ComplexDouble k2 ///Coefficient 2
								 ) const;///Computes the integral \int V(x) cos(k1*x) cos(k2*x)

  ComplexDouble FastSinIntegrate(ComplexDouble k1, ///Coefficient 1
								 ComplexDouble k2 ///Coefficient 2
								 ) const; ///Computes the integral \int V(x) sin(k1*x)sin(k2*x)

  void AddValue(Interval toAdd
				); ///Adds a value.

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
					  ) const; ///Integrate the basis functions over the potential to create V(k1, k2). 

  
  void Clear(); ///Clear the points in the potential.
  
  unsigned long GetPrecision(); ///Returns the precision.

  void SetPrecision(unsigned long value ///The value of the new precision.
					); ///Set the precision to the value given by 'value'.

private:
  void InitializeDefault(); ///Todo: replace with reading from file etc.
  list<Interval> PotentialPoints;
  double minX, maxX;
  unsigned long precision;
};



#endif
