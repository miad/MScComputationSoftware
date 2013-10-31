#include "LegendreRuleTest.hh"

bool LegendreRuleTest::TestCase1() const
{
  ///1 point
  vector<pair<double, double> > myVec = LegendreRule::GetRule(1);
  if(myVec.size() != 1 )
	return false;
  if( !DBL_EQUAL(myVec.front().first, 0) || ! DBL_EQUAL(myVec.front().second,2) )
	return false;
  return true;
}

bool LegendreRuleTest::TestCase2() const
{
  ///2 points
  vector<pair<double, double> > myVec = LegendreRule::GetRule(2);
  if(myVec.size() != 2 )
	return false;
  sort(myVec.begin(), myVec.end());

  if( !DBL_EQUAL(myVec[0].first, -sqrt(3.)/3.) && !DBL_EQUAL(myVec[0].second, 1. ))
	return false;
  if( !DBL_EQUAL(myVec[1].first, sqrt(3.)/3.) && !DBL_EQUAL(myVec[1].second, 1. ))
	return false;
  return true;
}


bool LegendreRuleTest::TestCase3() const
{
  ///3 points
  vector<pair<double, double> > myVec = LegendreRule::GetRule(3);
  if(myVec.size() != 3 )
	return false;
  sort(myVec.begin(), myVec.end());

  if( !DBL_EQUAL(myVec[0].first, -sqrt(3./5.)) && !DBL_EQUAL(myVec[0].second, 5./9. ))
	return false;
  if( !DBL_EQUAL(myVec[1].first, 0) && !DBL_EQUAL(myVec[1].second, 8./9. ))
	return false;
  if( !DBL_EQUAL(myVec[2].first, sqrt(3./5.)) && !DBL_EQUAL(myVec[2].second, 5./9. ))
	return false;

  return true;
}


bool LegendreRuleTest::TestCase4() const
{
  ///4 points
  vector<pair<double, double> > myVec = LegendreRule::GetRule(4);
  if(myVec.size() != 4 )
	return false;
  sort(myVec.begin(), myVec.end());

  if( !DBL_EQUAL(myVec[0].first, -sqrt((3.+2.*sqrt(6./5.))/7.)) && !DBL_EQUAL(myVec[0].second, (18.-sqrt(30.))/36.) )
	return false;

  if( !DBL_EQUAL(myVec[1].first, -sqrt((3.-2.*sqrt(6./5.))/7.)) && !DBL_EQUAL(myVec[1].second, (18.+sqrt(30.))/36.) )
	return false;

  if( !DBL_EQUAL(myVec[2].first, sqrt((3.-2.*sqrt(6./5.))/7.)) && !DBL_EQUAL(myVec[2].second, (18.+sqrt(30.))/36.) )
	return false;

  if( !DBL_EQUAL(myVec[3].first, sqrt((3.+2.*sqrt(6./5.))/7.)) && !DBL_EQUAL(myVec[3].second, (18.-sqrt(30.))/36.) )
	return false;

  return true;
}


int LegendreRuleTest::runUnitTests() const
{
  cout << "Running unit tests on LegendreRule...";
  cout << flush;
  if(!TestCase1())
    return 1;
  if(!TestCase2())
    return 2;
  if(!TestCase3())
    return 3;
  if(!TestCase4())
    return 4;

  cout << "done" << endl;
  return 0;
}
