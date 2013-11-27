#include "IntervalTest.hh"

int IntervalTest::TestCase1() const
{
  Interval myInterval(-1,2,5);
  if(myInterval.x1 != -1 || myInterval.x2 != 2 || myInterval.y != 5)
	return 1;
  Interval myInterval2(1,-1,4);
  if(myInterval2.x1 != -1 || myInterval2.x2 != 1 || myInterval2.y != 4)
	return 2;
  if(!myInterval.Overlaps(myInterval2))
	return 3;
  if(myInterval.Overlaps(Interval(2,3,4)))
	return 4;
  return 0;
}

int IntervalTest::TestCase2() const
{
  Interval I1 (1, 2, 3);
  Interval I2 (2, 3, 1);
  Interval I3 (1, 2, 3);
  Interval I4 (1, 2, 4);
  if( I1 > I2 )
	return 1;
  if( I1 == I2 )
	return 2;
  if(I1 != I3 )
	return 3;
  if( I1 == I4)
	return 4;
  if(I4 > I2 )
	return 5;
  if(I1 < I3 )
	return 6;
  return 0;
}

int IntervalTest::runUnitTests() const
{
  cout << "Running unit tests on Interval...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  int code2 = TestCase2();
  if(code2)
	return code2 + 100;
  cout << "done" << endl;
  return 0;
}
