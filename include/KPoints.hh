#ifndef KPoints_hh
#define KPoints_hh 1

#include <vector>
#include "Globals.hpp"
#include "RLException.hh"
#include <iostream>

using namespace std;

class KPoints
{
public:
  KPoints(unsigned int _nPoints, 
		  double _kCutoff, 
		  double _kMid, 
		  double _kDepth
		  ); ///Constructor.
  ~KPoints(); ///Destructor.
  vector< ComplexDouble > * GetPoints() const;
  unsigned int GetNPoints() const; ///Getter
  double GetKCutoff() const; ///Getter 
  double GetKMid() const; ///Getter 
  double GetKDepth() const; ///Getter
  ComplexDouble GetPoint(unsigned int i) const; /// Returns a specific point, with index i.
  ComplexDouble GetDeltaK(unsigned int i) const; /// Returns the local spacing at the point with index i.
  ComplexDouble GetStencilDeltaK() const;
private:
  vector< ComplexDouble > * kPoints;
  void ValidateArguments(unsigned int nPoints); /// Throws an exception if the arguments are not conforming with the expected requirements put on them. 
  
  double kCutoff; ///Maximum k-value before cutoff.
  double kMid; ///Mid k value for the Berggren curve.
  double kDepth; ///Depth for the Berggren curve.
};


#endif
