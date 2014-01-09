#include "LapackeEigenvalueSolverTest.hh"

int LapackeEigenvalueSolverTest::TestCase1() const
{
  ///Verify that the solver throws on a non-square matrix.
  CMatrix myMatrix(4,3);
  myMatrix.InitializeAll(ComplexDouble(4, 2));
  try
	{
	  EigenInformation myInfo = LapackeEigenvalueSolver::Solve(&myMatrix);
	}
  catch(RLException & ex)
	{
	  return 0;
	}
  return 1;
}

int LapackeEigenvalueSolverTest::TestCase2() const
{
  CMatrix myMatrix(3,3);
  myMatrix.InitializeAll(ComplexDouble(0, 0));
  for(int i = 0; i<3; ++i)
	myMatrix.Element(i, i) = ComplexDouble(i+1, i);
  EigenInformation myInfo = LapackeEigenvalueSolver::Solve(&myMatrix);
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

int LapackeEigenvalueSolverTest::TestCase3() const
{
  CMatrix myMatrix(3, 3);
  myMatrix.InitializeAll(ComplexDouble(0, 0));

  myMatrix.Element(0, 0) = ComplexDouble(0.814723686393179, 0.964888535199277);
  myMatrix.Element(0, 1) = ComplexDouble(0.913375856139019, 0.957166948242946);
  myMatrix.Element(0, 2) = ComplexDouble(0.278498218867048, 0.141886338627215);
  myMatrix.Element(1, 0) = ComplexDouble(0.905791937075619, 0.157613081677548);
  myMatrix.Element(1, 1) = ComplexDouble(0.632359246225410, 0.485375648722841);
  myMatrix.Element(1, 2) = ComplexDouble(0.546881519204984, 0.421761282626275);
  myMatrix.Element(2, 0) = ComplexDouble(0.126986816293506, 0.970592781760616);
  myMatrix.Element(2, 1) = ComplexDouble(0.097540404999410, 0.800280468888800);
  myMatrix.Element(2, 2) = ComplexDouble(0.957506835434298, 0.915735525189067);

  EigenInformation myInfo = LapackeEigenvalueSolver::Solve(&myMatrix);


  /*
  cout << endl;  
  for(int i = 0; i<3; ++i)
	{
	  printf("λ = %+1.4lf%+1.4lf → V = [", real(myInfo.Eigenvalues[i]), imag(myInfo.Eigenvalues[i]));
	  for(int j = 0; j<3; ++j)
		{
		  printf("%+1.4lf%+1.4lf", real(myInfo.Eigenvectors[i][j]), imag(myInfo.Eigenvectors[i][j]));
		  if(j<2)
			printf(", ");
		}
	  printf("].\n");
	}
  */

  if(myInfo.Eigenvalues.size() != 3)
	return 1;
  if(myInfo.Eigenvectors.size() != 3)
	return 2;
  for(int i = 0; i<3; ++i)
	{
	  if(DBL_EQUAL(myInfo.Eigenvalues[i],  ComplexDouble(1.879212794739092, 1.903170649787818)))
		{
		  if(myInfo.Eigenvectors[i].size() != 3)
			return 3;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][0], ComplexDouble(0.485282228366207, -0.296708052870155)))
			return 4;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][1], ComplexDouble(0.338848858887965, - 0.308475666349372)))
			{
			  return 5;
			}
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][2], ComplexDouble(0.683000515588799, 0.000000000000000)))
			return 6;
		  
		}
	  else if(DBL_EQUAL(myInfo.Eigenvalues[i], ComplexDouble(-0.293596624026318, 0.186556816736693)))
		{
		  if(myInfo.Eigenvectors[i].size() != 3)
			return 7;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][0], ComplexDouble(-0.673526180545294, - 0.162803324616888)))
			return 8;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][1], ComplexDouble(0.710590551892770, 0.000000000000000)))
			return 9;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][2], ComplexDouble(-0.047837784181806, + 0.112384053935867)))
			return 10;
		}
	  else if(DBL_EQUAL(myInfo.Eigenvalues[i], ComplexDouble(0.818973597340112, + 0.276272242586674)))
		{
		  if(myInfo.Eigenvectors[i].size() != 3)
			return 11;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][0], ComplexDouble(-0.464235641088859, - 0.209300657923136)))
			return 12;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][1], ComplexDouble(-0.080089459524766, + 0.308785039594677)))
			return 13;
		  if(!EIGEN_EQUAL(myInfo.Eigenvectors[i][2], ComplexDouble(0.799322201575374, + 0.000000000000000)))
			return 14;
		}
	  else
		{
		  return 40;
		}

	}
  return 0;
}


int LapackeEigenvalueSolverTest::runUnitTests() const
{
  cout << "Running unit tests on LapackeEigenvalueSolver...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;
  code1 = TestCase2();
  if(code1)
	return 100+code1;
  code1 = TestCase3();
  if(code1)
	return 200 + code1;
  cout << "done" << endl;
  return 0;
}


bool LapackeEigenvalueSolverTest::CplexCompare(const ComplexDouble & c1, const ComplexDouble & c2)
{
  if(real(c1) > real(c2))
	return true;
  else if(real(c2) > real(c1))
	return false;
  else if(imag(c1) > imag(c2))
	return true;
  return false;
}
