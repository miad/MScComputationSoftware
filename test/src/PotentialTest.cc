#include "PotentialTest.hh"

int PotentialTest::TestCase1() const
{
  Potential myPotential;
  myPotential.AddValue(-1,0,1);
  myPotential.AddValue(2,3,3);
  if(myPotential.Evaluate(-2) != 0)
	return 1;
  if(myPotential.Evaluate(-0.5) != 1)
	return 2;
  if(myPotential.Evaluate(1.5) != 0)
	return 3;
  if(myPotential.Evaluate(2.5) != 3)
	return 4;
  if(!DBL_EQUAL(myPotential.FastExpIntegrate(ComplexDouble(0,0)),ComplexDouble(4, 0)))
	return 5;
  if(!DBL_EQUAL(myPotential.FastExpIntegrate(ComplexDouble(PI,0)),ComplexDouble(0, 4/PI)))
	return 6;
  try ///Should trow if we try to add an overlap.
	{
	  myPotential.AddValue(-.5, 1.5, 1);
	  return 7;
	}
  catch(RLException & ex)
	{
	  ///nothing here.
	}
  return 0;
}

int PotentialTest::runUnitTests() const
{
  cout << "Running unit tests on Potential...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  cout << "done" << endl;
  return 0;
}
