#ifndef ParametrizedCurve_hh
#define ParametrizedCurve_hh 1

#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include "Globals.hpp"
#include "LegendreRule.hh"


using namespace std;

/** 
	A parametrized curve in the complex plane, consisting of straight lines between points in the complex plane.
	Possible generalization: allow a curve other than a straight line (given by an arbitrary function?) between any two points. Priority: low.
	The ParametrizedCurve class can also associate, and keep track of, Gauss-Legendre quadrature rules for each segment. 

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

  void AddGLPoints(unsigned int numberOfPoints ///Number of GL points on the curve.
				   );

  ComplexDouble Evaluate(double parameterValue
						 ) const; ///Evaluate the curve at a specific position, given by parameterValue. Complexity O(log(n)), where n is the number of points on the curve.

  ComplexDouble SegmentEvaluate(unsigned int segment, 
								double parameterValue
								) const; ///Evaluate the specified segment for the specified parameterValue. O(1). Note that parameterValue goes from start to stop over the segment.

  double GetStart() const; ///Returns the start value of the parameter for the curve.
  double GetStop() const; ///Returns the stop value of the parameter for the curve.
  unsigned int GetNumberOfValues() const; ///Returns the number of values added to the curve.
  unsigned int GetNumberOfSegments() const; ///Returns the number of segments (the number of values minus one).

  ComplexDouble GetSegmentDerivative(unsigned int segment ///Choose which segment vector to return.
								  ) const; ///Returns the difference between the endpoints of the curve, divided by (stop-start) to normalize. The derivative (which is constant) of the segment line.
  
  double GetLength() const; ///Returns the curve length.
  void ComputeParameterValues(); ///Returns the length of the curve. Guaranteed not to be negative. Complexity O(n), where n is the number of parameters added to the curve. Can probably be improved to O(1) using some cleaverness.


  const vector<pair<ComplexDouble, ComplexDouble> > * GetSegmentRule(unsigned int segment ///Segment number.
													  ) const; ///Returns the GL rule for a specific segment. Implicitly calls ComputeGaussLegendre if this has never been done before.

  void ComputeGaussLegendre(); ///Computes GL rules for all segments. Throws an exception if all segments has not been associated with a rule number.

  void Clear(); ///Clears the object of all information.

  unsigned int GetTotalNumberOfGLPoints(); ///Returns the total number of GL points.

  unsigned int SegmentIndexFromGLNumber(unsigned int val
										); ///Returns the segment index from a GL number.


  const pair<ComplexDouble, ComplexDouble> * GetRulePoint(unsigned int segment, 
											 unsigned int GLpoint
											 ) const;

  ComplexDouble GetRuleValue(unsigned int segment, 
							 unsigned int GLpoint
							 ) const;

  ComplexDouble GetRuleWeight(unsigned int segment, 
					   unsigned int GLpoint
					   ) const;



private:

  unsigned int totalNumberOfGLPoints;

  vector<pair<double, ComplexDouble> > ParametrizedCurvePoints; ///The points on the curve, together with the parameter value to which the points correspond. The parameter value is rescaled every time a value is added.

  vector<unsigned int> numberOfGLPoints; ///Number of Gauss-Legendre points, used when constructing the GL rules.

  vector< vector<pair<ComplexDouble, ComplexDouble > > > gaussLegendreValues; ///The outer vector contains a vector for each segment. Each segment vector contains a number of Gauss-Legendre quadrature points (the value on the curve as well as the weight), that is, the curve value on these points.


  static bool ComparePairs(const pair<double, ComplexDouble> & p1, ///First pair
						   const pair<double, ComplexDouble> & p2 ///Second pair
						   );///Comparision function.

  double start, stop; ///The positions on the curve for start and stop.
  double length; ///The curve length.
};


#endif
