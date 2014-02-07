#include "ArpackEigenvalueSolver.hh"


//ArpackError::ErrorCode ArpackError::code = NO_ERRORS;


EigenInformation * ArpackEigenvalueSolver::Solve(CMatrix * toSolve, uint numberOfEigenvalues, ComplexDouble shift)
{

  if( ! toSolve->IsSquare() )
	{
	  throw RLException("Can only find eigenvalues of square matrices.");
	}
  
  /*
  uint n = toSolve->Rows();
  ComplexDouble * a = toSolve->GetArray();
  ARComplexDouble * matr = new ARComplexDouble[n*n];
  for(uint i = 0; i<n; ++i)
	matr[i] = a[i];
  


  ARdsNonSymMatrix<ARComplexDouble, ARfloat> A(n*n, a);

  ARComplexDouble sigma = shift;

  
  ARluCompStdEig<ARfloat> Prob(numberOfEigenvalues, A, sigma);
  
  Prob.FindEigenvectors();

  if(Prob.ConvergedEigenvalues() != (int)numberOfEigenvalues)
	{
	  throw RLException("Only %d eigenvalues converged, %d were requested.\n",Prob.ConvergedEigenvalues(), numberOfEigenvalues);
	}
  
  EigenInformation * toReturn = new EigenInformation();
  toReturn->Eigenvalues.resize(Prob.ConvergedEigenvalues());
  toReturn->Eigenvectors.resize(Prob.ConvergedEigenvalues(), vector<ComplexDouble>(n, 0.0) );
  for(int i = 0; i<Prob.ConvergedEigenvalues(); ++i)
	{
	  toReturn->Eigenvalues.at(i) = Prob.Eigenvalue(i);
	  for(uint j = 0; j<n; ++j)
		{
		  toReturn->Eigenvectors.at(i).at(j) = Prob.Eigenvector(i, j);
		}
	}

  delete[] matr; matr = NULL;

  return toReturn;
  */
  return NULL;
}
