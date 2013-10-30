#include "ParametrizedCurve.hh"

ParametrizedCurve::ParametrizedCurve(double _start, double _stop)
  :start(_start), stop(_stop)
{
  if(stop < start)
	swap(start, stop);
}

ParametrizedCurve::~ParametrizedCurve()
{

}

void ParametrizedCurve::AddValue(ComplexDouble val, unsigned int pos)
{
  length = -1;
  list<pair<double, ComplexDouble> >::iterator posToInsert = ParametrizedCurvePoints.begin(); 
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
  
  for(list<pair<double, ComplexDouble> >::iterator it = ParametrizedCurvePoints.begin(); it!=ParametrizedCurvePoints.end(); ++it)
	{
	  length += abs(lastValue-it->second);
	  it->first = length;
	  lastValue = it->second;
	}
		
  for(list<pair<double, ComplexDouble> >::iterator it = ParametrizedCurvePoints.begin(); it!=ParametrizedCurvePoints.end(); ++it)
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

  list<pair<double, ComplexDouble> >::const_iterator it = upper_bound(ParametrizedCurvePoints.begin(), ParametrizedCurvePoints.end(), make_pair(x, ComplexDouble(0)), ComparePairs);
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

ComplexDouble ParametrizedCurve::GetStart() const
{
  return start;
}

ComplexDouble ParametrizedCurve::GetStop() const
	{
  return stop;
}

unsigned int ParametrizedCurve::GetNumberOfValues() const
	{
  return ParametrizedCurvePoints.size();
}

bool ParametrizedCurve::ComparePairs(const pair<double, ComplexDouble> & p1, const pair<double, ComplexDouble> & p2)
{
  return p1.first < p2.first;
}
