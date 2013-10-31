#ifndef Potential_hh
#define Potential_hh 1

#include <list>
#include <set>
#include "Interval.hh"
#include "RLException.hh"
#include "Globals.hpp"

#ifndef EPS
#define EPS (1E-9)
#endif


using namespace std;

/** 
	Stepwise potential
 */
class Potential
{
public:
  Potential(); ///Constructor. \todo Implement reading potential from file.
  ~Potential(); ///Destructor.
  void AddValue(double x1, ///The x1 value to add.
				double x2, ///The x2 value to add.
				double y ///The y value to add.
				); ///Adds a value to the stepwise potential.

  ComplexDouble FastExpIntegrate(ComplexDouble exponentVal ///The exponent in the fast integration.
								 ); ///Fast integration of e^{i*exponentVal*x} over the potential. Returns the integral over the entire potential.

  void AddValue(Interval toAdd
				); ///Adds a value.

  double Evaluate(double x ///The point to evaluate the potential value at.
				  ) const; ///Evaluate a specific point in this potential.
  double GetMinX() const; ///Returns the minimum x value for which the potential is nonzero.
  double GetMaxX() const; ///Returns the maximum x value for which the potential is nonzero.
  unsigned int GetNumberOfValues() const;
private:
  void InitializeDefault(); ///Todo: replace with reading from file etc.
  list<Interval> PotentialPoints;
  double minX, maxX;
};

double PotentialValue(double x);

void InitializePotential();

void ClearPotential();

void SetPoint(int n, float x1, float x2, float y);



#endif
