#include "HarmonicBasisFunctionTest.hh"

int HarmonicBasisFunctionTest::TestCase1() const
{
  //Using Legendre rule.

  HermiteEvaluator::Init(1, true);
  double hbar = 7.63823302; ///Magic value. Unit kB*µK*µs
  double mass = 723.453025853; ///Total mass in special mass units.
  double omega = 4.270E-3;
  double xmin = 4;
  SpecificUnits myUnits(hbar, mass, "","",1.0);
  
  ParametrizedPotential myPotential("x", vector<pair<string, double> >(), -5, 5);

  vector<Potential*> potV; potV.push_back(&myPotential);

  HarmonicBasisFunction basis(xmin, omega, &potV, &myUnits);

  vector<pair<double, double> > rule = LegendreRule::GetRule(300, -20, 25);

  for(uint n = 0; n<NMAX_CERT_ONE; ++n)
	{
	  for(uint m = 0; m<NMAX_CERT_ONE; ++m)
		{
		  long double sum = 0;
		  for(vector<pair<double, double> >::const_iterator it = rule.begin(); it!=rule.end(); ++it)
			{
			  sum += it->second * basis.Eval(n, it->first) * basis.Eval(m, it->first);
			}
		  if(!GL_EQ(KRONECKER_DELTA(n, m), sum) )
			{
			  return 1+n+m;
			}
		}
	}
  
  return 0;
}


int HarmonicBasisFunctionTest::TestCase2() const
{
  //Using Hermite rule.

  HermiteEvaluator::Init(1, true);
  double hbar = 7.63823302; ///Magic value. Unit kB*µK*µs
  double mass = 723.453025853; ///Total mass in special mass units.
  double omega = 4.270E-3;
  double xmin = 4;
  SpecificUnits myUnits(hbar, mass, "µm","µK*kB", 0.0);

  vector<pair<string, double> > args;
  args.push_back(make_pair("mass", mass));
  args.push_back(make_pair("omega", omega));
  args.push_back(make_pair("xmin", xmin));
  
  
  ParametrizedPotential myPotential("0.5*mass*pow(omega, 2)*pow(x-xmin,2)", args, -50000, 50000);
  vector<Potential*> potV; potV.push_back(&myPotential);


  HarmonicBasisFunction basis(xmin, omega, &potV, &myUnits);

  vector<pair<double, double> > rule = HermiteRule::GetRule(300, xmin, mass*omega/hbar); //myUnits.GetMassOverLambda2()*omega/(myUnits.GetHbarTimesLambda()));
  for(uint n = 0; n<NMAX_CERT_TWO; ++n)
	{
	  for(uint m = 0; m<NMAX_CERT_TWO; ++m)
		{
		  double sum = 0;
		  
		  for(vector<pair<double, double> >::const_iterator it = rule.begin(); it!=rule.end(); ++it)
			{
			  sum += it->second * basis.EvalNonExponentPart(n, it->first) * basis.EvalNonExponentPart(m, it->first);
			}
		  if( !GL_EQ(KRONECKER_DELTA(n, m), sum) )
			{
			  cout << "n " << n << " m " << m << " sum " << sum << endl;
			  return 1+n+m;
			}
		}
	}

  return 0;
}



int HarmonicBasisFunctionTest::TestCase3() const
{
  HermiteEvaluator::Init(1, true);
  double hbar = 7.63823302; ///Magic value. Unit kB*µK*µs
  double mass = 723.453025853; ///Total mass in special mass units.
  double omega = 4.270E-3;
  double xmin = 4;
  SpecificUnits myUnits(hbar, mass, "","",1.0);

  vector<pair<string, double> > args;
  args.push_back(make_pair("mass", mass));
  args.push_back(make_pair("omega", omega));
  args.push_back(make_pair("xmin", xmin));
  
  
  ParametrizedPotential myPotential("0.5*mass*pow(omega, 2.0)*pow((x-xmin),2.0)", args, -50000, 50000);
  vector<Potential*> potV; potV.push_back(&myPotential);


  HarmonicBasisFunction basis(xmin, omega, &potV, &myUnits);

  vector<pair<double, double> > rule = HermiteRule::GetRule(300, xmin, mass*omega/hbar);
  for(uint n = 0; n<NMAX_CERT; ++n)
	{
	  for(uint m = 0; m<NMAX_CERT; ++m)
		{
		  double sum = 0;
		  
		  for(vector<pair<double, double> >::const_iterator it = rule.begin(); it!=rule.end(); ++it)
			{
			  sum += it->second * basis.EvalNonExponentPart(n, it->first) * basis.EvalNonExponentPart(m, it->first) * myPotential.Evaluate(it->first);
			}
		  ///Compare with the matrix elements for X^2 in the harmonic oscillator case.
		  double expect = 
			hbar*omega/(4.0)* (
								  KRONECKER_DELTA(n, m) * (2*m+1) + 
								  KRONECKER_DELTA(m, n-2)* sqrt(n*(n-1)) + 
								  KRONECKER_DELTA(m, n+2) * sqrt((n+1)*(n+2))
								  );
		  if(!GL_EQ(expect, sum) )
			return 1+n+m;
		}
	}

  return 0;
}



int HarmonicBasisFunctionTest::TestCase4() const
{
  HermiteEvaluator::Init(1, true);
  double hbar = 7.63823302; ///Magic value. Unit kB*µK*µs
  double mass = 723.453025853; ///Total mass in special mass units.
  double omega = 4.270E-3;
  double xmin = 4;
  SpecificUnits myUnits(hbar, mass, "","",1.0);

  vector<pair<string, double> > args;
  args.push_back(make_pair("mass", mass));
  args.push_back(make_pair("omega", omega));
  args.push_back(make_pair("xmin", xmin));
  
  
  ParametrizedPotential myPotential("0.5*mass*pow(omega, 2.0)*pow((x-xmin),2.0)-1.0", args, -50000, 50000);
  vector<Potential*> potV; potV.push_back(&myPotential);


  HarmonicBasisFunction basis(xmin, omega, &potV, &myUnits);

  vector<pair<double, double> > rule = HermiteRule::GetRule(300, xmin, mass*omega/hbar);
  for(uint n = 0; n<NMAX_CERT; ++n)
	{
	  for(uint m = 0; m<NMAX_CERT; ++m)
		{
		  if(basis.DiffIntegrate(n, m) != 0)
			return n + m + 1;
		}
	}

  return 0;
}

int HarmonicBasisFunctionTest::TestCase5() const
{
  HermiteEvaluator::Init(1, true);
  double hbar = 7.63823302; ///Magic value. Unit kB*µK*µs
  double mass = 723.453025853; ///Total mass in special mass units.
  double omega = 4.270E-3;
  double xmin = 4;
  SpecificUnits myUnits(hbar, mass, "","",1.0);

  vector<pair<string, double> > args;
  args.push_back(make_pair("mass", mass));
  args.push_back(make_pair("omega", omega));
  args.push_back(make_pair("xmin", xmin));
  
  
  ParametrizedPotential myPotential("0.5*mass*pow(omega, 2.0)*pow((x-xmin),2.0) - 1.0", args, -50000, 50000);
  vector<Potential*> potV; potV.push_back(&myPotential);


  HarmonicBasisFunction basis(xmin, omega, &potV, &myUnits);

  vector<pair<double, double> > rule = HermiteRule::GetRule(300, xmin, mass*omega/hbar);
  for(uint n = 0; n<NMAX_CERT; ++n)
	{
	  if(! DBL_EQUAL(basis.GetEigenEnergy(n), (0.5L+n)*hbar*omega - 1.0))
		return n;
	}

  return 0;
}



int HarmonicBasisFunctionTest::runUnitTests() const
{
  cout << "Running unit tests on HarmonicBasisFunction...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return 1;
  int code2 = TestCase2();
  if(code2)
	return 2;

  int code3 = TestCase3(); 
  if(code3)
	return 3;

  int code4 = TestCase4();
  if(code4)
	return 4;

  int code5 = TestCase5();
  if(code5)
	return 5;

  cout << "done" << endl;
  return 0;
}
