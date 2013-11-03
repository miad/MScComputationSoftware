#include "EigenvalueSolverTest.hh"

int EigenvalueSolverTest::TestCase1() const
{
  ///Verify that the solver throws on a non-square matrix.
  CMatrix myMatrix(4,3);
  myMatrix.InitializeAll(ComplexDouble(4, 2));
  try
	{
	  EigenInformation myInfo = EigenvalueSolver::Solve(&myMatrix);
	}
  catch(RLException & ex)
	{
	  return 0;
	}
  return 1;
}

int EigenvalueSolverTest::TestCase2() const
{
  CMatrix myMatrix(3,3);
  myMatrix.InitializeAll(ComplexDouble(0, 0));
  for(int i = 0; i<3; ++i)
	myMatrix.Element(i, i) = ComplexDouble(i+1, i);
  EigenInformation myInfo = EigenvalueSolver::Solve(&myMatrix);
  sort(myInfo.Eigenvalues.rbegin(), myInfo.Eigenvalues.rend(), CplexCompare);

  for(int i = 0; i<3; ++i)
	{
	  if(!DBL_EQUAL(myInfo.Eigenvalues[i], ComplexDouble(i+1, i)))
		return i+1;
	}

  for(int j = 0; j<3; ++j)
	{
	  for(int i = 0; i<3; ++i)
		{
		  if(!DBL_EQUAL(myInfo.Eigenvectors[i][j], ComplexDouble((i==j),0) ))
			return 10+10*i+j;
		}
	}
  return 0;
}


int EigenvalueSolverTest::runUnitTests() const
{
  cout << "Running unit tests on EigenvalueSolver...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  code1 = TestCase2();
  if(code1)
	return 100+code1;
  cout << "done" << endl;
  return 0;
}


bool EigenvalueSolverTest::CplexCompare(const ComplexDouble & c1, const ComplexDouble & c2)
{
  if(real(c1) > real(c2))
	return true;
  else if(real(c2) > real(c1))
	return false;
  else if(imag(c1) > imag(c2))
	return true;
  return false;
}
