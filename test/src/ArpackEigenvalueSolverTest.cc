#include "ArpackEigenvalueSolverTest.hh"

int ArpackEigenvalueSolverTest::TestCase1() const
{
  ///Verify that the solver throws on a non-square matrix.
  CMatrix myMatrix(4,3);
  myMatrix.InitializeAll(ComplexDouble(4, 2));
  try
	{
	  EigenInformation * myInfo = ArpackEigenvalueSolver::Solve(&myMatrix, 0, 0.0, true);
	  delete myInfo;
	}
  catch(RLException & ex)
	{
	  return 0;
	}
  return 1;
}

int ArpackEigenvalueSolverTest::TestCase2() const
{
  srand(time(NULL));
  CMatrix myMatrix(1000,1000);
  for(int i = 0; i<1000; ++i)
	for(int j = 0; j<1000; ++j)
	  {
		myMatrix.Element(i, j) = ComplexDouble(rand() % 100, rand() % 100);
	  }

  EigenInformation * myInfo = ArpackEigenvalueSolver::Solve(&myMatrix, 3, 0.0, true);

  delete myInfo;
  return 0;
}

int ArpackEigenvalueSolverTest::runUnitTests() const
{
  cout << "Running unit tests on ArpackEigenvalueSolver...";
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


bool ArpackEigenvalueSolverTest::CplexCompare(const ComplexDouble & c1, const ComplexDouble & c2)
{
  if(real(c1) > real(c2))
	return true;
  else if(real(c2) > real(c1))
	return false;
  else if(imag(c1) > imag(c2))
	return true;
  return false;
}
