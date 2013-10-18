#include "rfacLinearInterpolationTest.hh"

bool rfacLinearInterpolationTest::rfacLinearInterpolation_ConstructsCorrectly_AssertTrue() const
{
  vector<double > coordinates; coordinates.push_back(0); coordinates.push_back(1);
  vector<double> values; values.push_back(1); values.push_back(2);
  rfacLinearInterpolation myInterpolation(coordinates, values);

  vector<pair<double, double> > data;
  data.push_back(make_pair(1,2)); data.push_back(make_pair(3,4));
  rfacLinearInterpolation myInterpolation2(data);

  return true;
}


bool rfacLinearInterpolationTest::rfacLinearInterpolation_ThrowsOnTooSmallInput_AssertTrue() const
{
  int thrown = 0;
  try
    {
      vector<double > coordinates; coordinates.push_back(0);
      vector<double> values; values.push_back(1);
      rfacLinearInterpolation myInterpolation(coordinates, values);
    }
  catch(rfacException e)
    {
      ++thrown;
    }
  assert(thrown==1);

  try
    {
      vector<pair<double, double> > data;
      data.push_back(make_pair(1,2));
      rfacLinearInterpolation myInterpolation2(data);
    }
  catch(rfacException e)
    {
      ++thrown;
    }
  assert(thrown==2);
  return thrown==2;
}

rfacLinearInterpolation rfacLinearInterpolationTest::initializeDefault1() const
{
  /*
    The governing equations are as following:
    range -3 to 0 y=-0.33333x +1
    range 0 to 1 y=1x +1
    range 1 to 3 y=1.5x +0.5
    range 3 to 4 y=8x +-19
    range 4 to 9.5 y=-1.9273x +20.7090909091
  */
  vector<pair<double, double> > data;
  data.push_back(make_pair(-3,2));
  data.push_back(make_pair(0,1));
  data.push_back(make_pair(1,2));
  data.push_back(make_pair(3,5));
  data.push_back(make_pair(4,13));
  data.push_back(make_pair(9.5,2.4));
  rfacLinearInterpolation toReturn(data);
  return toReturn;
}

bool rfacLinearInterpolationTest::rfacLinearInterpolation_PerformsCorrectFitI_AssertTrue(double i) const
{
  if(i<=0)
    {
      assert(abs(initializeDefault1().Eval(i)-(-0.333333333333*i+1))<EPSILON);
    }
  if(i>=0 && i<=1)
    {
      assert(abs(initializeDefault1().Eval(i)-(i+1))<EPSILON);
    }
  if(i>=1 && i<=3)
    {
      assert(abs(initializeDefault1().Eval(i)-(1.5*i+0.5))<EPSILON);
    }
  if(i>=3 && i<=4)
    {
      assert(abs(initializeDefault1().Eval(i)-(8*i-19))<EPSILON);
    }
  if(i>=4)
    {
      assert(abs(initializeDefault1().Eval(i)-(-1.92727272727272727272727*i+20.709090909090909090))<EPSILON);
    }


  return true;
}

int rfacLinearInterpolationTest::runUnitTests() const
{
  cout << "Running unit tests on rfacLinearInterpolation...";
  cout << flush;
  if(!rfacLinearInterpolation_ConstructsCorrectly_AssertTrue())
    return 1;
  if(!rfacLinearInterpolation_ThrowsOnTooSmallInput_AssertTrue())
    return 2;
  for(double i = -10; i<=100; i+=0.01)
    {
      if(!rfacLinearInterpolation_PerformsCorrectFitI_AssertTrue(i))
	return (int)i*10;
    }
  cout << "done" << endl;
  return 0;
}
