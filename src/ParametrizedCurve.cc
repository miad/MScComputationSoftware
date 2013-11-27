#include "ParametrizedCurve.hh"

ParametrizedCurve::ParametrizedCurve(double _start, double _stop)
  :start(_start), stop(_stop)
{
  if(stop < start)
	swap(start, stop);
  if(stop == start)
	throw RLException("Interval can not have zero length.");
  totalNumberOfGLPoints = 0;
}

ParametrizedCurve::~ParametrizedCurve()
{

}

unsigned int ParametrizedCurve::GetTotalNumberOfGLPoints()
{
  return totalNumberOfGLPoints;
}

unsigned int ParametrizedCurve::SegmentIndexFromGLNumber(unsigned int val)
{
  if(val >= totalNumberOfGLPoints)
	{
	  throw RLException("Invalid GL number to get segment index from.");
	}

  for(unsigned int i = 0; i<numberOfGLPoints.size(); ++i)
	{
	  if(cumNumberOfGLPoints[i] > val)
		{
		  return i;
		}
	}
  throw RLException("Could not identify index from GL number. This is an internal inconsistency.");
}

void ParametrizedCurve::AddValue(ComplexDouble val, unsigned int pos)
{
  length = -1;
  vector<pair<double, ComplexDouble> >::iterator posToInsert = ParametrizedCurvePoints.begin(); 
  advance(posToInsert, pos);
  ParametrizedCurvePoints.insert(posToInsert, make_pair(0.0, val));

  ComputeParameterValues();
}

void ParametrizedCurve::AddGLPoints(unsigned int numberOfPoints)
{
  if(numberOfGLPoints.size() > GetNumberOfSegments())
	{
	  throw RLException("Tried to add more GL points than there are segments.");
	}

  numberOfGLPoints.push_back(numberOfPoints);
  totalNumberOfGLPoints += numberOfPoints;
  cumNumberOfGLPoints.push_back(totalNumberOfGLPoints);
}

void ParametrizedCurve::ComputeGaussLegendre()
{
  if(numberOfGLPoints.size() != GetNumberOfSegments())
	{
	  throw RLException("Cannot compute Gauss-Legendre points: number of GL points added does not correspond to the number of segments added. %d segments but %d values.", GetNumberOfSegments(),  numberOfGLPoints.size());
	}
  gaussLegendreValues.clear();
  int segment = 0;
  for(vector<unsigned int>::const_iterator it = numberOfGLPoints.begin(); it!=numberOfGLPoints.end(); ++it)
	{
	  gaussLegendreValues.push_back(vector<pair<ComplexDouble, ComplexDouble> >());
	  vector<pair<double, double> > rule = LegendreRule::GetRule(*it, start, stop);
	  for(vector<pair<double, double> >::const_iterator ip = rule.begin(); ip != rule.end(); ++ip)
		{
		  gaussLegendreValues.back().push_back(make_pair(SegmentEvaluate(segment, ip->first),
														 ip->second*GetSegmentDerivative(segment)));
		}
	  ++segment;
	}
}

void ParametrizedCurve::Clear()
{
  totalNumberOfGLPoints = 0;
  ParametrizedCurvePoints.clear();
  numberOfGLPoints.clear();
  cumNumberOfGLPoints.clear();
  gaussLegendreValues.clear();
  length = -1;
}

void ParametrizedCurve::AddValue(ComplexDouble val)
{
  length = -1;
  ParametrizedCurvePoints.push_back(make_pair(0, val));
  ComputeParameterValues();
}

void ParametrizedCurve::ComputeParameterValues()
{
  length = 0;
  if(ParametrizedCurvePoints.empty())
	{
	  return;
	}

  ComplexDouble lastValue = ParametrizedCurvePoints.front().second;
  
  for(vector<pair<double, ComplexDouble> >::iterator it = ParametrizedCurvePoints.begin(); it!=ParametrizedCurvePoints.end(); ++it)
	{
	  length += abs(lastValue-it->second);
	  it->first = length;
	  lastValue = it->second;
	}
		
  for(vector<pair<double, ComplexDouble> >::iterator it = ParametrizedCurvePoints.begin(); it!=ParametrizedCurvePoints.end(); ++it)
	{
	  it->first = it->first/length*(stop-start) + start;
	}
}

double ParametrizedCurve::GetLength() const
{
  if( length < 0 )
	{
	  throw RLException("Length was less than 0.");
	}
  return length;
}

ComplexDouble ParametrizedCurve::Evaluate(double x) const
	{
  if(x < start || x > stop)
	throw RLException("Parameter for evaluation was out of acceptable range: parameter was %lf, acceptable range [%lf, %lf].", x, start, stop);

  vector<pair<double, ComplexDouble> >::const_iterator it = upper_bound(ParametrizedCurvePoints.begin(), ParametrizedCurvePoints.end(), make_pair(x, ComplexDouble(0)), ComparePairs);
  if(it == ParametrizedCurvePoints.begin())
	{
	  throw RLException("Error: comparision fail in parametrized curve. This should not happen.");
	}
  ComplexDouble p1val = it->second;
  double p1param = it->first;
  --it;
  ComplexDouble p0val = it->second;
  double p0param = it->first;

  //One-point formula.
  return p0val + (p1val-p0val)/(p1param-p0param)*(x-p0param);
}

ComplexDouble ParametrizedCurve::SegmentEvaluate(unsigned int segment, double parameterValue) const
{
  if( segment + 1 >= ParametrizedCurvePoints.size())
	throw RLException("Cannot evaluate segment %d: out of bounds.", segment);
  ComplexDouble p1val = ParametrizedCurvePoints[segment+1].second;
  ComplexDouble p0val = ParametrizedCurvePoints[segment].second;
  return p0val + (p1val-p0val)/ComplexDouble((stop-start),0)*(parameterValue-start);
}

ComplexDouble ParametrizedCurve::GetSegmentDerivative(unsigned int segment) const
{
  if( segment + 1 >= ParametrizedCurvePoints.size() )
	throw RLException("Cannot get segment vector for segment %d: out of bounds.", segment);
  return (ParametrizedCurvePoints[segment+1].second - ParametrizedCurvePoints[segment].second)/(ComplexDouble(stop-start,0));
}

const vector<pair<ComplexDouble, ComplexDouble> > * ParametrizedCurve::GetSegmentRule(unsigned int segment) const
{
  if(segment >= gaussLegendreValues.size())
	{
	  throw RLException("Asked for GL rule for segment %d, but there are only %d rules computed. Did you forget a call to ComputeGaussLegendre() ?\n", segment, gaussLegendreValues.size());
	}
  return &gaussLegendreValues[segment];
}

const pair<ComplexDouble, ComplexDouble> * ParametrizedCurve::GetRulePoint(unsigned int segment, unsigned int GLpoint) const
{
  unsigned int offset = (segment > 0)?cumNumberOfGLPoints[segment-1]:0;
  if(offset > GLpoint)
	{
	  throw RLException("Internal inconsistency: offset > GLPoint.");
	}
  if( segment > GetNumberOfSegments() )
	{
	  throw RLException("Invalid segment number: requested %d but there are %d segments.", segment, GetNumberOfSegments());
	}
  if( GLpoint - offset >= gaussLegendreValues[segment].size() )
	{
	  throw RLException("Invalid GL point number on segment %d: requested %d but max is %d.", segment, GLpoint, gaussLegendreValues[segment].size());
	}

  return &gaussLegendreValues[segment][GLpoint-offset];
}

ComplexDouble ParametrizedCurve::GetRuleValue(unsigned int segment, unsigned int GLpoint) const
{
  return GetRulePoint(segment, GLpoint)->first;
}

ComplexDouble ParametrizedCurve::GetRuleWeight(unsigned int segment, unsigned int GLpoint) const
{
  return GetRulePoint(segment, GLpoint)->second;
}

double ParametrizedCurve::GetStart() const
{
  return start;
}

double ParametrizedCurve::GetStop() const
	{
  return stop;
}

unsigned int ParametrizedCurve::GetNumberOfValues() const
{
  return ParametrizedCurvePoints.size();
}

unsigned int ParametrizedCurve::GetNumberOfSegments() const
{
  return ParametrizedCurvePoints.size() - 1;
}

bool ParametrizedCurve::ComparePairs(const pair<double, ComplexDouble> & p1, const pair<double, ComplexDouble> & p2)
{
  return p1.first < p2.first;
}
