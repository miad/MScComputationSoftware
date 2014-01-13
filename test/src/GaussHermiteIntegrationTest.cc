#include "GaussHermiteIntegrationTest.hh"

int GaussHermiteIntegrationTest::TestCase1() const
{
  for(int a = -5; a<5; ++a)
	{
	  for(uint b = 1; b<10; ++b)
		{
		  vector<pair<double, double> > myRule = HermiteRule::GetRule(200, a, b);
		  double myVal0 = 0;
		  for(vector<pair<double, double> >::const_iterator it = myRule.begin(); it!=myRule.end(); ++it)
			{
			  myVal0 += it->second ;
			}
		  if(!DBL_EQUAL(myVal0, sqrt(PI/b)))
			return 1;
		}
	}
  return 0;
}


int GaussHermiteIntegrationTest::TestCase2() const
{
  double xmin = 4;
  vector<pair<double, double> > rule = HermiteRule::GetRule(200, xmin, PI);
  double sum = 0;
  for(vector<pair<double, double> >::const_iterator it = rule.begin(); it!=rule.end(); ++it)
	{
	  sum += it->second;
	}
  if(!DBL_EQUAL(sum, 1))
	return 1;

  return 0;

}

int GaussHermiteIntegrationTest::runUnitTests() const
{
  cout << "Running integration tests on Gauss-Hermite rules...";
  cout << flush;
  int retcode = TestCase1();
  if(retcode)
	return retcode;

  retcode = TestCase1();
  if(retcode)
	return 1+retcode;


  retcode = TestCase2();
  if(retcode)
	return 2+retcode;


  cout << "done" << endl;
  return 0;
}
