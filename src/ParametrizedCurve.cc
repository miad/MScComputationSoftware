#include "ParametrizedCurve.hh"

ParametrizedCurve::ParametrizedCurve(double _start, double _stop)
  :start(_start), stop(_stop)
{
  if(stop < start)
	swap(start, stop);
  if(stop == start)
	throw RLException("Interval can not have zero length.");
}

ParametrizedCurve::~ParametrizedCurve()
{

}

void ParametrizedCurve::AddValue(ComplexDouble val, unsigned int pos)
{
  length = -1;
  vector<pair<double, ComplexDouble> >::iterator posToInsert = ParametrizedCurvePoints.begin(); 
  advance(posToInsert, pos);
  ParametrizedCurvePoints.insert(posToInsert, make_pair(0.0, val));
  ComputeParameterValues();
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
