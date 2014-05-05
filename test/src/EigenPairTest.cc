#include "EigenPairTest.hh"

int EigenPairTest::TestCase1() const
{
  EigenPair p1, p2;
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(6, 10);
  if(!(p1 < p2))
	return 1;
  if(p2 < p1)
	return 2;
  return 0;
}

int EigenPairTest::TestCase2() const
{
  EigenPair p1, p2;
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 5);
  if(!(p1 < p2))
	return 1;
  if(p2 < p1)
	return 2;
  return 0;

}

int EigenPairTest::TestCase3() const
{
  EigenPair p1, p2;
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  if((p1 < p2))
	return 1;
  if((p2 < p1))
	return 2;
  return 0;

}

int EigenPairTest::TestCase4() const
{
  EigenPair p1(2), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  if(!(p1 < p2))
	return 1;
  if((p2 < p1))
	return 2;
  return 0;

}

int EigenPairTest::TestCase5() const
{
  EigenPair p1(5), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  if((p1 < p2))
	return 1;
  if(!(p2 < p1))
	return 2;
  return 0;

}



int EigenPairTest::TestCase6() const
{
  EigenPair p1(3), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  if((p1 < p2))
	return 1;
  if((p2 < p1))
	return 2;
  return 0;

}

int EigenPairTest::TestCase7() const
{
  EigenPair p1(3), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  p1.Eigenvector[0] = ComplexDouble(1,1);
  p1.Eigenvector[1] = ComplexDouble(1,2);
  p2.Eigenvector[0] = ComplexDouble(1,1);
  p2.Eigenvector[1] = ComplexDouble(1,3);

  if(!(p1 < p2))
	return 1;
  if((p2 < p1))
	return 2;
  return 0;
}

int EigenPairTest::TestCase8() const
{
  EigenPair p1(3), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  p1.Eigenvector[0] = ComplexDouble(1,1);
  p1.Eigenvector[1] = ComplexDouble(1,4);
  p2.Eigenvector[0] = ComplexDouble(1,1);
  p2.Eigenvector[1] = ComplexDouble(1,3);

  if((p1 < p2))
	return 1;
  if(!(p2 < p1))
	return 2;
  return 0;
}

int EigenPairTest::TestCase9() const
{
  EigenPair p1(3), p2(3);
  p1.Eigenvalue = ComplexDouble(3, 4);
  p2.Eigenvalue = ComplexDouble(3, 4);
  p1.Eigenvector[0] = ComplexDouble(1,1);
  p1.Eigenvector[1] = ComplexDouble(1,4);
  p2.Eigenvector[0] = ComplexDouble(1,1);
  p2.Eigenvector[1] = ComplexDouble(1,4);

  if((p1 < p2))
	return 1;
  if((p2 < p1))
	return 2;
  return 0;
}




int EigenPairTest::runUnitTests() const
{
  cout << "Running unit tests on EigenPair...";
  cout << flush;
  if(TestCase1())
	return 1;
  if(TestCase2())
	return 2;
  if(TestCase3())
	return 3;
  if(TestCase4())
	return 4;
  if(TestCase5())
	return 5;
  if(TestCase6())
	return 6;
  if(TestCase7())
	return 7;
  if(TestCase8())
	return 8;
  if(TestCase9())
	return 9;

  cout << "done" << endl;
  return 0;
}
