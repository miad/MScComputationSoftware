#ifndef ParametrizedCurve_hh
#define ParametrizedCurve_hh 1

#include <list>
#include <set>
#include <algorithm>
#include <utility>
#include "Globals.hpp"



using namespace std;

/** 
	A parametrized curve in the complex plane, consisting of straight lines between points in the complex plane.
	Possible generalization: allow a curve other than a straight line (given by an arbitrary function?) between any two points. Priority: low.

 */
class ParametrizedCurve
{
public:
  ParametrizedCurve(double _start = -1, ///Initial value of the parameter.
					double _stop = 1 ///Final value of the parameter. Must not be == _start. If _stop < _start, these values will be swapped.
					); ///Constructor. 

  ~ParametrizedCurve(); ///Destructor.

  void AddValue(ComplexDouble val, ///Value to insert
				unsigned int pos ///Position to insert to.
				); ///Adds a value at a certain position.

  void AddValue(ComplexDouble val ///Value to add.
				) ;///Adds a value to the end of the curve. Complexity O(n) seen from caller..


  ComplexDouble Evaluate(double parameterValue
						 ) const; ///Evaluate the curve at a specific position, given by parameterValue. Complexity O(log(n)), where n is the number of points on the curve.

  ComplexDouble GetStart() const; ///Returns the start value of the parameter for the curve.
  ComplexDouble GetStop() const; ///Returns the stop value of the parameter for the curve.
  unsigned int GetNumberOfValues() const; ///Returns the number of values added to the curve.
  double GetLength() const; ///Returns the curve length.
  void ComputeParameterValues(); ///Returns the length of the curve. Guaranteed not to be negative. Complexity O(n), where n is the number of parameters added to the curve. Can probably be improved to O(1) using some cleaverness.

private:
  list<pair<double, ComplexDouble> > ParametrizedCurvePoints; ///The points on the curve, together with the parameter value to which the points correspond. The parameter value is rescaled every time a value is added.

  static bool ComparePairs(const pair<double, ComplexDouble> & p1, ///First pair
						   const pair<double, ComplexDouble> & p2 ///Second pair
						   );///Comparision function.

  double start, stop; ///The positions on the curve for start and stop.
  double length; ///The curve length.
};


#endif
