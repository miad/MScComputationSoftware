#ifndef Potential_hh
#define Potential_hh 1

///Abstract class describing any 1D potential.
class Potential
{
public:
  ///Abstract functions. Any potential is supposed to implement these.
  double Evaluate(double x ///The point to evaluate the potential value at.
				  ) const = 0; ///Evaluate a specific point in this potential.
  double GetMinX() const = 0; ///Returns the minimum x value for which the potential is nonzero.
  double GetMaxX() const = 0; ///Returns the maximum x value for which the potential is nonzero.
  unsigned long GetPrecision() const = 0; ///Returns the precision with which evaluations are performed. The implementation of this is arbitrary and potential-dependent, but a suggestion is some relation to the number of Gauss-Legendre points used when integrating.
  void SetPrecision(unsigned long value ///The value to set the precision to.
					); ///Sets the precision for the integration to Value.
  ComplexDouble BasisIntegrate(BasisFunction & b1, ///First basis function.
					  BasisFunction & b2, ///Second basis function.
					  ComplexDouble & k1, ///First k-value.
					  ComplexDouble & k2 ///Second k-value.
					  ) const = 0; ///Integrate the basis functions over the potential to create V(k1, k2). NOTE: This is supposed to be a FAST implementation, and faster than a naive one based on Evaluate since the Potential class can optimize itself dependent on how it is implemented. 
}

#endif
