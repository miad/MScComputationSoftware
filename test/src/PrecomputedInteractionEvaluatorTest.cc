#include "PrecomputedInteractionEvaluatorTest.hh"

int PrecomputedInteractionEvaluatorTest::TestCase1() const
{
  PrecomputedInteractionEvaluatorExposer myExposer;
  if(myExposer.Permutations_X(4, 4, 4, 4) != 1)
	return 1;
  if(myExposer.Permutations_X(4, 4, 4, 3) != 4)
	return 2;
  if(myExposer.Permutations_X(4, 4, 3, 4) != 4)
	return 3;
  if(myExposer.Permutations_X(4, 3, 4, 4) != 4)
	return 4;
  if(myExposer.Permutations_X(3, 4, 4, 4) != 4)
	return 5;
  if(myExposer.Permutations_X(1, 1, 4, 4) != 6)
	return 6;
  if(myExposer.Permutations_X(1, 4, 1, 4) != 6)
	return 7;
  if(myExposer.Permutations_X(1, 4, 4, 1) != 6)
	return 8;
  if(myExposer.Permutations_X(1, 2, 3, 4) != 24)
	return 9;
  if(myExposer.Permutations_X(1, 1, 3, 4) != 12)
	return 10;
  if(myExposer.Permutations_X(1, 2, 2, 4) != 12)
	return 11;
  if(myExposer.Permutations_X(1, 2, 3, 3) != 12)
	return 12;
  if(myExposer.Permutations_X(1, 2, 3, 1) != 12)
	return 13;

  return 0;
}

int PrecomputedInteractionEvaluatorTest::TestCase2() const
{
  PrecomputedInteractionEvaluatorExposer myExposer;
  uint nmax = 7;
  vector<vector<vector<ComplexDouble> > > PsiB(2, vector<vector<ComplexDouble> >(10, vector<ComplexDouble>(nmax, 1.0) ) );
  vector<vector<vector<vector<double> > > > Vnnnn(nmax, vector<vector<vector<double> > >(nmax, vector<vector<double> >(nmax, vector<double>(nmax, 1.0) ) ));
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 3, 4, Vnnnn, PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 3, 4, Vnnnn, PsiB, nmax)
				) 
	 )
	return 1;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 3, 3, Vnnnn, PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 3, 3, Vnnnn, PsiB, nmax)
				) 
	 )
	return 2;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 2, 1, 2, Vnnnn, PsiB, nmax), 
				myExposer.ComputeElement_R(1, 2, 1, 2, Vnnnn, PsiB, nmax)
				) 
	 )
	return 3;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(1, 1, 1, 4, Vnnnn, PsiB, nmax), 
				myExposer.ComputeElement_R(1, 1, 1, 4, Vnnnn, PsiB, nmax)
				) 
	 )
	return 4;
  if(!DBL_EQUAL(myExposer.ComputeElement_X(2, 2, 2, 2, Vnnnn, PsiB, nmax), 
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
  int code1 = TestCase1();
  if(code1)
    return code1;
  int code2 = TestCase2();
  if(code2)
	return 100+code2;
  cout << "done" << endl;
  return 0;
}

