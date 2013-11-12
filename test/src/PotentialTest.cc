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

int PotentialTest::TestCase2() const
{
  Potential myPotential("testsample/potentialtest1.tsv");
  list<Interval> myIs = myPotential.GetPotentialPoints();
  int ptr = 0;
  for(list<Interval>::iterator it = myIs.begin(); it!=myIs.end(); ++it)
	{
	  switch(ptr)
		{
		case 0:
		  if(!DBL_EQUAL(it->x1,1) || !DBL_EQUAL(it->x2,3) || !DBL_EQUAL(it->y,5))
			return 1;
		  break;
		case 1:
		  if(!DBL_EQUAL(it->x1,5) || !DBL_EQUAL(it->x2,6) || !DBL_EQUAL(it->y,2))
			return 2;
		  break;
		case 2:
		  if(!DBL_EQUAL(it->x1,-1) || !DBL_EQUAL(it->x2,1) || !DBL_EQUAL(it->y,0.5))
			return 3;
		  break;
		default:
		  return 4;
		}
	  ++ptr;
	}
  return 0;
}


int PotentialTest::TestCase3() const
{
  Potential myPotential;
  myPotential.AddValue(-1,0,1);
  myPotential.AddValue(2,3,3);
  
  if(!DBL_EQUAL(myPotential.FastCosIntegrate(0, 0),4.))
	return 1;
  if(!DBL_EQUAL(myPotential.FastSinIntegrate(0, 0),0.))
	return 2;
  if(!DBL_EQUAL(myPotential.FastCosIntegrate(1., 2.),-0.36224364268112103165))
	return 3;
  if(!DBL_EQUAL(myPotential.FastCosIntegrate(1.4, -2.12),-0.6069969616186261))
	return 4;
  if(!DBL_EQUAL(myPotential.FastCosIntegrate(1, 1),2.5853646045381722078))
	return 5;
  if(!DBL_EQUAL(myPotential.FastSinIntegrate(1.4, -2.12),-.8558303084692681))
	return 6;
  if(!DBL_EQUAL(myPotential.FastSinIntegrate(1,1),1.4146353954618277922))
	return 7;
  if(!DBL_EQUAL(myPotential.FastCosIntegrate(1.5, -1.5),2.369286993063653))
	return 8;
  if(!DBL_EQUAL(myPotential.FastSinIntegrate(1.4,-1.4),-1.144181287850408))
	return 9;


  return 0;
}


int PotentialTest::runUnitTests() const
{
  cout << "Running unit tests on Potential...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  code1 = TestCase2();
  if(code1)
	return 10+code1;
  code1 = TestCase3();
  if(code1)
	return 100 + code1;


  cout << "done" << endl;
  return 0;
}
