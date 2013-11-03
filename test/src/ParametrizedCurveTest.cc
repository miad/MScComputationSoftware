#include "ParametrizedCurveTest.hh"

bool ParametrizedCurveTest::TestCase1() const
{
  ///Sets the boundaries correctly.
  ParametrizedCurve myCurve(-4.3,3.22);
  
  if(! DBL_EQUAL(myCurve.GetStart(), -4.3))
	return false;
  if(! DBL_EQUAL(myCurve.GetStop(), 3.22))
	return false;
  return true;
}

bool ParametrizedCurveTest::TestCase2() const
{
  ///Reverses boundaries if they are wrongly directed.
  ParametrizedCurve myCurve(4.3,-3.22);
  
  if(! DBL_EQUAL(myCurve.GetStart(), -3.22))
	return false;
  if(! DBL_EQUAL(myCurve.GetStop(), 4.3))
	return false;
  return true;
}

bool ParametrizedCurveTest::TestCase3() const
{
  ///Parametrizes a trivial line correctly: k(t) = t, -10 < t < 10
  ParametrizedCurve myCurve(-10, 10);
  
  myCurve.AddValue(-10);
  myCurve.AddValue(10);
  
  for(double i = -10; i<=10; i+=0.1)
	if(! DBL_EQUAL(myCurve.Evaluate(i), ComplexDouble(i)))
	  return false;
  return true;
}

bool ParametrizedCurveTest::TestCase4() const
{
  ///Parametrizes a trivial line correctly: k(t) = -t, -10 < t < 10
  ParametrizedCurve myCurve(-10, 10);
  
  myCurve.AddValue(10);
  myCurve.AddValue(-10);

  
  for(double i = -10; i<=10; i+=0.1)
	{
	  if(! DBL_EQUAL(myCurve.Evaluate(i), ComplexDouble(-i)))
	  return false;
	}
  return true;
}

bool ParametrizedCurveTest::TestCase5() const
{
  ///Parametrizes a trivial line correctly: k(t) = t, -10 < t < 10
  ParametrizedCurve myCurve(-10, 10);

  myCurve.AddValue(-10, 0);  
  myCurve.AddValue(10, 1);
  myCurve.AddValue(5, 1);

  
  for(double i = -10; i<=10; i+=0.1)
	{
	  if(! DBL_EQUAL(myCurve.Evaluate(i), ComplexDouble(i,0)))
		return false;
	}
  return true;
}


int ParametrizedCurveTest::TestCase6() const
{
  //Full-fledged test. 
  ParametrizedCurve myCurve(0, 28);

  myCurve.AddValue(ComplexDouble(-10, 1));  
  myCurve.AddValue(ComplexDouble(-6, 1));
  myCurve.AddValue(ComplexDouble(-3, 5));
  myCurve.AddValue(ComplexDouble(3, -3));
  myCurve.AddValue(ComplexDouble(7, 0));
  myCurve.AddValue(ComplexDouble(11, 0));

  if(myCurve.GetNumberOfValues() != 6)
	return 1;

  if(! DBL_EQUAL(myCurve.GetLength(), ComplexDouble(28, 0)))
	return 2;

  if(! DBL_EQUAL(myCurve.Evaluate(0), ComplexDouble(-10, 1)))
	return 3;
  if(! DBL_EQUAL(myCurve.Evaluate(1), ComplexDouble(-9, 1)))
	return 4;
  if(! DBL_EQUAL(myCurve.Evaluate(2), ComplexDouble(-8, 1)))
	return 5;
  if(! DBL_EQUAL(myCurve.Evaluate(3), ComplexDouble(-7, 1)))
	return 6;
  if(! DBL_EQUAL(myCurve.Evaluate(4), ComplexDouble(-6, 1)))
	return 7;
  if(! DBL_EQUAL(myCurve.Evaluate(5), ComplexDouble(-6+3./5., 1+4./5.)))
	return 8;
  if(! DBL_EQUAL(myCurve.Evaluate(6), ComplexDouble(-6+2*3./5., 1+2*4./5.)))
	return 9;
  if(! DBL_EQUAL(myCurve.Evaluate(7), ComplexDouble(-6+3*3./5., 1+3*4./5.)))
	return 10;
  if(! DBL_EQUAL(myCurve.Evaluate(8), ComplexDouble(-6+4*3./5., 1+4*4./5.)))
	return 11;
  if(! DBL_EQUAL(myCurve.Evaluate(9), ComplexDouble(-3, 5)))
	return 12;
  if(! DBL_EQUAL(myCurve.Evaluate(10), ComplexDouble(-3+3./5., 5-4./5.)))
	return 13;
  if(! DBL_EQUAL(myCurve.Evaluate(11), ComplexDouble(-3+2*3./5., 5-2*4./5.)))
	return 14;
  if(! DBL_EQUAL(myCurve.Evaluate(12), ComplexDouble(-3+3*3./5., 5-3*4./5.)))
	return 15;
  if(! DBL_EQUAL(myCurve.Evaluate(13), ComplexDouble(-3+4*3./5., 5-4*4./5.)))
	return 16;
  if(! DBL_EQUAL(myCurve.Evaluate(14), ComplexDouble(-3+5*3./5., 5-5*4./5.)))
	return 17;
  if(! DBL_EQUAL(myCurve.Evaluate(15), ComplexDouble(-3+6*3./5., 5-6*4./5.)))
	return 18;
  if(! DBL_EQUAL(myCurve.Evaluate(16), ComplexDouble(-3+7*3./5., 5-7*4./5.)))
	return 19;
  if(! DBL_EQUAL(myCurve.Evaluate(17), ComplexDouble(-3+8*3./5., 5-8*4./5.)))
	return 20;
  if(! DBL_EQUAL(myCurve.Evaluate(18), ComplexDouble(-3+9*3./5., 5-9*4./5.)))
	return 21;
  if(! DBL_EQUAL(myCurve.Evaluate(19), ComplexDouble(3, -3)))
	return 22;
  if(! DBL_EQUAL(myCurve.Evaluate(20), ComplexDouble(3 + 4./5., -3 + 3./5.)))
	return 23;
  if(! DBL_EQUAL(myCurve.Evaluate(21), ComplexDouble(3 + 2*4./5., -3 + 2*3./5.)))
	return 24;
  if(! DBL_EQUAL(myCurve.Evaluate(22), ComplexDouble(3 + 3*4./5., -3 + 3*3./5.)))
	return 25;
  if(! DBL_EQUAL(myCurve.Evaluate(23), ComplexDouble(3 + 4*4./5., -3 + 4*3./5.)))
	return 26;
  if(! DBL_EQUAL(myCurve.Evaluate(24), ComplexDouble(7, 0)))
	return 27;
  if(! DBL_EQUAL(myCurve.Evaluate(25), ComplexDouble(8, 0)))
	return 28;
  if(! DBL_EQUAL(myCurve.Evaluate(26), ComplexDouble(9, 0)))
	return 29;
  if(! DBL_EQUAL(myCurve.Evaluate(27), ComplexDouble(10, 0)))
	return 30;
  if(! DBL_EQUAL(myCurve.Evaluate(28), ComplexDouble(11, 0)))
	return 31;

  return 0;
}


int ParametrizedCurveTest::TestCase7() const
{
  //Full-fledged test. 
  ParametrizedCurve myCurve(-2,2);

  myCurve.AddValue(ComplexDouble(-10, 1));  
  myCurve.AddValue(ComplexDouble(-6, 1));
  myCurve.AddValue(ComplexDouble(-3, 5));
  myCurve.AddValue(ComplexDouble(3, -3));
  myCurve.AddValue(ComplexDouble(7, 0));
  myCurve.AddValue(ComplexDouble(11, 0));

  if(myCurve.GetNumberOfSegments() != 5)
	return 1;

  if(! DBL_EQUAL(myCurve.GetSegmentDerivative(0), ComplexDouble(1., 0.)))
	return 2;
  if(! DBL_EQUAL(myCurve.GetSegmentDerivative(1), ComplexDouble(3./4., 4./4.)))
	return 3;
  if(! DBL_EQUAL(myCurve.GetSegmentDerivative(2), ComplexDouble(6./4., -8./4.)))
	return 4;
  if(! DBL_EQUAL(myCurve.GetSegmentDerivative(3), ComplexDouble(4./4., 3./4.)))
	return 5;
  if(! DBL_EQUAL(myCurve.GetSegmentDerivative(4), ComplexDouble(4./4., 0.)))
	return 6;
  try
	{
	  ComplexDouble vask = myCurve.GetSegmentDerivative(5);
	  (void)vask;
	  return 7;
	}
  catch(RLException & ex){}





  if(! DBL_EQUAL(myCurve.SegmentEvaluate(0,-2), ComplexDouble(-10, 1)))
	return 8;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(0,-1), ComplexDouble(-9, 1)))
	return 9;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(0,0), ComplexDouble(-8, 1)))
	return 10;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(0,1), ComplexDouble(-7, 1)))
	return 11;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(0,2), ComplexDouble(-6, 1)))
	return 12;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(1,-2+4./5.), ComplexDouble(-6+3./5., 1+4./5.)))
	return 13;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(1,-2+2*4./5.), ComplexDouble(-6+2*3./5., 1+2*4./5.)))
	return 14;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(1,-2+3*4./5.), ComplexDouble(-6+3*3./5., 1+3*4./5.)))
	return 15;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(1,-2+4*4./5.), ComplexDouble(-6+4*3./5., 1+4*4./5.)))
	return 16;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(1,2), ComplexDouble(-3, 5)))
	return 17;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+1*4./10.), ComplexDouble(-3+3./5., 5-4./5.)))
	return 18;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+2*4./10.), ComplexDouble(-3+2*3./5., 5-2*4./5.)))
	return 19;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+3*4./10.), ComplexDouble(-3+3*3./5., 5-3*4./5.)))
	return 20;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+4*4./10.), ComplexDouble(-3+4*3./5., 5-4*4./5.)))
	return 21;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+5*4./10.), ComplexDouble(-3+5*3./5., 5-5*4./5.)))
	return 22;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+6*4./10.), ComplexDouble(-3+6*3./5., 5-6*4./5.)))
	return 23;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+7*4./10.), ComplexDouble(-3+7*3./5., 5-7*4./5.)))
	return 24;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+8*4./10.), ComplexDouble(-3+8*3./5., 5-8*4./5.)))
	return 25;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+9*4./10.), ComplexDouble(-3+9*3./5., 5-9*4./5.)))
	return 26;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(2,-2+10*4./10.), ComplexDouble(3, -3)))
	return 27;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(3, -2+1*4./5.), ComplexDouble(3 + 4./5., -3 + 3./5.)))
	return 28;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(3, -2+2*4./5.), ComplexDouble(3 + 2*4./5., -3 + 2*3./5.)))
	return 29;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(3, -2+3*4./5.), ComplexDouble(3 + 3*4./5., -3 + 3*3./5.)))
	return 30;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(3, -2+4*4./5.), ComplexDouble(3 + 4*4./5., -3 + 4*3./5.)))
	return 31;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(3, -2+5*4./5.), ComplexDouble(7, 0)))
	return 32;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(4, -1), ComplexDouble(8, 0)))
	return 33;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(4, 0), ComplexDouble(9, 0)))
	return 34;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(4, 1), ComplexDouble(10, 0)))
	return 35;
  if(! DBL_EQUAL(myCurve.SegmentEvaluate(4, 2), ComplexDouble(11, 0)))
	return 36;

  return 0;
}





int ParametrizedCurveTest::runUnitTests() const
{
  cout << "Running unit tests on ParametrizedCurve...";
  cout << flush;
  if(!TestCase1())
    return 1;
  if(!TestCase2())
    return 2;
  if(!TestCase3())
    return 3;
  if(!TestCase4())
    return 4;
  if(!TestCase5())
    return 5;

  int retCode = TestCase6();
  if(retCode)
    return (10+retCode);
  retCode = TestCase7();
  if(retCode)
	return 100 + retCode;

  cout << "done" << endl;
  return 0;
}
