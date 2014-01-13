#ifndef Potential_hh
#define Potential_hh 1

#include "Globals.hpp"
#include "BasisFunction.hh"
#include <utility>
#include <vector>


enum PotentialType
  {
	TypePiecewiseConstant = 1,
	TypeParametrized = 0
  };



///Abstract class describing any 1D potential.
class Potential
{
public:

  virtual ~Potential()  = 0 ; ///Virtual destructor.
  
  ///Abstract functions. Any potential is supposed to implement these.
  virtual double Evaluate(double x ///The point to evaluate the potential value at.
				  ) const = 0; ///Evaluate a specific point in this potential.

  virtual double GetMinX() const = 0; ///Returns the minimum x value for which the potential is nonzero.
  virtual double GetMaxX() const = 0; ///Returns the maximum x value for which the potential is nonzero.
  virtual unsigned long GetPrecision() const = 0; ///Returns the precision with which evaluations are performed. The implementation of this is arbitrary and potential-dependent, but a suggestion is some relation to the number of Gauss-Legendre points used when integrating.
  virtual void SetPrecision(unsigned long value ///The value to set the precision to.
							) = 0; ///Sets the precision for the integration to Value.
  virtual ComplexDouble BasisIntegrate(BasisFunction & b1, ///First basis function.
									   BasisFunction & b2, ///Second basis function.
									   ComplexDouble & k1, ///First k-value.
									   ComplexDouble & k2 ///Second k-value.
									   ) = 0; ///Integrate the basis functions over the potential to create V(k1, k2). NOTE: This is supposed to be a FAST implementation, and faster than a naive one based on Evaluate since the Potential class can optimize itself dependent on how it is implemented. 

  virtual vector<pair<double, double> > GetPlottingPoints() const = 0; ///Returns a vector with points (implemented as pairs) which may be used for plotting. The points should be such that if they are connected by lines in order, they will give a realistic picture of the potential.

  virtual vector<pair<double, double> > GetPrecisionPoints() const = 0; ///Returns points suitable for plotting indicating the precision.

  virtual PotentialType GetType() const = 0;

};

#endif
