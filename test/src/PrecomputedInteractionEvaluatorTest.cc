#include "PrecomputedInteractionEvaluatorTest.hh"



int PrecomputedInteractionEvaluatorTest::TestCase1() const
{
  PrecomputedInteractionEvaluatorExposer myExposer;
  uint nmax = 7;
  vector<vector<vector<ComplexDouble> > > PsiB(2, vector<vector<ComplexDouble> >(10, vector<ComplexDouble>(nmax, 1.0) ) );
  vector<vector<vector<vector<double> > > > Vnnnn(nmax, vector<vector<vector<double> > >(nmax, vector<vector<double> >(nmax, vector<double>(nmax, 1.0) ) ));
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 3, 4, &Vnnnn, &PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 3, 4, Vnnnn, PsiB, nmax)
				) 
	 )
	return 1;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 3, 3, &Vnnnn, &PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 3, 3, Vnnnn, PsiB, nmax)
				) 
	 )
	return 2;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 1, 2, &Vnnnn, &PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 1, 2, Vnnnn, PsiB, nmax)
				) 
	 )
	return 3;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 1, 1, 4, &Vnnnn, &PsiB, nmax), 
				myExposer.ComputeElement_R(1, 1, 1, 4, Vnnnn, PsiB, nmax)
				) 
	 )
	return 4;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(2, 2, 2, 2, &Vnnnn, &PsiB, nmax), 
				myExposer.ComputeElement_R(2, 2, 2, 2, Vnnnn, PsiB, nmax)
				) 
	 )
	return 5;
  return 0;
}

int PrecomputedInteractionEvaluatorTest::runUnitTests() const
{
  cout << "Running unit tests on PrecomputedInteractionEvaluator...";
  cout << flush;
  int code2 = TestCase1();
  if(code2)
	return 100+code2;
  cout << "done" << endl;
  return 0;
}

