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

int IntervalTest::runUnitTests() const
{
  cout << "Running unit tests on Interval...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  cout << "done" << endl;
  return 0;
}
