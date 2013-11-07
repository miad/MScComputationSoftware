#include "GaussLegendreIntegrationTest.hh"

int GaussLegendreIntegrationTest::TestCase1() const
{
  ///Using the Cauchy theorem.
  vector<pair<double, double> > myVec = LegendreRule::GetRule(40);
  ParametrizedCurve kCurve;
  kCurve.AddValue(ComplexDouble(0.0,1.0));
  kCurve.AddValue(ComplexDouble(1.0,1.0));
  kCurve.AddValue(ComplexDouble(5.0,-2.0));
  kCurve.AddValue(ComplexDouble(6.0,4.2));
  kCurve.AddValue(ComplexDouble(6.4,5.0));
  kCurve.AddValue(ComplexDouble(-2.3,-3.3));
  kCurve.AddValue(ComplexDouble(0.0,1.0));
  
  for(int fun = 0; fun<4; ++fun)
	{
	  ComplexDouble integralValue = 0.;
	  for(unsigned int j = 0; j<kCurve.GetNumberOfSegments(); ++j)
		{
		  for(unsigned int i = 0; i<myVec.size(); ++i)
			{
			  
			  /*		  cout << "j= " << j << " seg-deriv: " << kCurve.GetSegmentDerivative(j) << endl;
						  cout << "eval x={0.3 0.6 0.9}: " << kCurve.SegmentEvaluate(j, 0.3) << " " << 
						  kCurve.SegmentEvaluate(j, 0.6) << " " << kCurve.SegmentEvaluate(j, 0.9) << endl;
			  */
			  ComplexDouble z = kCurve.SegmentEvaluate(j, myVec[i].first);
			  ComplexDouble f = 0;
			  switch(fun)
				{			  
				case 0: //constant.
				  f = 1;
				  break;
				case 1: //f(z) = z
				  f = z;
				  break;
				case 2: //f(z) = exp(z)
				  f = exp(z);
				  break;
				case 3:
				  f = cos(z) + ComplexDouble(3, 2) + exp(z)*3.;
				  break;
				}
			  integralValue += 
				kCurve.GetSegmentDerivative(j)*
				myVec[i].second*f;
			}
			  
		}
	  if(!DBL_EQUAL(integralValue, 0.))
		return fun + 1;
	}
	  
  return 0;
}

int GaussLegendreIntegrationTest::TestCase2() const
{
  ///Using the Cauchy theorem.
  vector<pair<double, double> > myVec = LegendreRule::GetRule(40);
  ParametrizedCurve kCurve;
  ///Random simple curve in complex plane.
  kCurve.AddValue(ComplexDouble(0.0,1.0));
  kCurve.AddValue(ComplexDouble(1.0,1.0));
  kCurve.AddValue(ComplexDouble(5.0,-2.0));
  kCurve.AddValue(ComplexDouble(6.0,4.2));
  kCurve.AddValue(ComplexDouble(6.4,5.0));
  kCurve.AddValue(ComplexDouble(7.0, -6.0));
  kCurve.AddValue(ComplexDouble(-2.3,-3.3));
  kCurve.AddValue(ComplexDouble(0.0,1.0));
  
  for(int fun = 0; fun<4; ++fun)
	{
	  ComplexDouble integralValue = 0.;
	  for(unsigned int j = 0; j<kCurve.GetNumberOfSegments(); ++j)
		{
		  for(unsigned int i = 0; i<myVec.size(); ++i)
			{
			  
			  /*		  cout << "j= " << j << " seg-deriv: " << kCurve.GetSegmentDerivative(j) << endl;
						  cout << "eval x={0.3 0.6 0.9}: " << kCurve.SegmentEvaluate(j, 0.3) << " " << 
						  kCurve.SegmentEvaluate(j, 0.6) << " " << kCurve.SegmentEvaluate(j, 0.9) << endl;
			  */
			  ComplexDouble z = kCurve.SegmentEvaluate(j, myVec[i].first);
			  ComplexDouble f = 0;
			  ComplexDouble point = ComplexDouble(1, -2);
			  switch(fun)
				{			  
				case 0: //constant.
				  f = 1.0/(z-point);
				  break;
				case 1: //f(z) = z
				  f = z/(z-point);
				  break;
				case 2: //f(z) = exp(z)
				  f = exp(z)/(z-point);
				  break;
				case 3:
				  f = (cos(z) + ComplexDouble(3, 2) + exp(z)*3.)/(z-point);
				  break;
				}
			  integralValue += 
				kCurve.GetSegmentDerivative(j)*
				myVec[i].second*f;
			}
			  
		}
	  switch(fun)
		{
		case 0:
		  if(!DBL_EQUAL(integralValue, -2*PI*ComplexDouble(0.0,1.0)))
			return fun + 1;
		  break;
		case 1:
		  if(!DBL_EQUAL(integralValue, -2*PI*ComplexDouble(0.0,1.0)*ComplexDouble(1, -2)))
			return fun + 1;
		  break;
		case 2:
		  if(!DBL_EQUAL(integralValue, -2*PI*ComplexDouble(0.0,1.0)*exp(ComplexDouble(1,-2))))
			return fun + 1;
		  break;
		case 3:
		  ComplexDouble z=ComplexDouble(1, -2);
		  if(!DBL_EQUAL(integralValue, -2*PI*ComplexDouble(0.0,1.0)*(cos(z) + ComplexDouble(3, 2) + exp(z)*3.)))
			return fun + 1;
		  break;
		}

	}
	  
  return 0;
}


int GaussLegendreIntegrationTest::runUnitTests() const
{
  cout << "Running integration tests on Gauss-Legendre Integration algorithm...";
  cout << flush;
  int retcode = TestCase1();
  if(retcode)
    return retcode;
  retcode = TestCase2();
  if(retcode)
	return 100 + retcode;

  cout << "done" << endl;
  return 0;
}
