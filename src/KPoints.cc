#include "KPoints.hh"


KPoints::KPoints(unsigned int nPoints, double _kCutoff, double _kMid, double _kDepth)
  : kCutoff(_kCutoff), kMid(_kMid), kDepth(_kDepth)
{
  kPoints = new vector< ComplexDouble >();
  double spacing = kCutoff / nPoints;
  
  foru(i, 2*(int)nPoints+1)
	{
	  double realPart = -1*kCutoff + spacing*i;
	  double imagPart = 0;
	  if( abs(realPart) > 0  && abs(realPart) <= kMid  && kMid > 0)
		{
		  realPart = SGN(realPart)*( -1* realPart * kDepth / kMid );
		} 
	  else if( abs(realPart) > kMid  && abs(realPart) <= 2*kMid  && kMid > 0)
		{
		  imagPart = SGN(realPart)*( -2*kDepth + realPart * kDepth / kMid );
		} 
	  kPoints->push_back(ComplexDouble(realPart, imagPart));
	}
}

void KPoints::ValidateArguments(unsigned int nPoints)
{
  if ( nPoints == 0)
	throw basisException("KPoints: Number of points asked for must be > 0.");
  if ( kMid < 0 )
	throw basisException("KPoints: kMid must be >= 0.");
  if ( kDepth < 0 )
	throw basisException("KPoints: kDepth must be >= 0.");
  if ( ( kMid > 0 && kDepth == 0 ) || (kMid == 0 && kDepth > 0 ) )
	throw basisException("KPoints: Either both of kMid and kDepth are 0, or none.");
  if ( kCutoff <= 2 * kMid )
	throw basisException("KPoints: kCutoff > 2 * kMid is required.");
}

KPoints::~KPoints()
{

}

vector<ComplexDouble> * KPoints::GetPoints() const
{
  return kPoints;
}

double KPoints::GetKCutoff() const
{
  return kCutoff;
}

double KPoints::GetKMid() const
{
  return kMid;
}

double KPoints::GetKDepth() const
{
  return kDepth;
}
