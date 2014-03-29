#include "ArpackEigenvalueSolver.hh"


ArpackError::ErrorCode ArpackError::code = NO_ERRORS;


EigenInformation * ArpackEigenvalueSolver::Solve(CMatrix * toSolve, ulong numberOfEigenvalues, ComplexDouble shift, bool findEigenvectors)
{


  if( ! toSolve->IsSquare() )
	{
	  throw RLException("Can only find eigenvalues of square matrices, dimensions %ld x %ld", toSolve->Rows(), toSolve->Columns());
	}

  if(numberOfEigenvalues == 0)
	return new EigenInformation();

  
  ulong n = toSolve->Rows();
  ComplexDouble * a = toSolve->GetArray();
  
  ARdsNonSymMatrix<ARComplexDouble, ARfloat> A(n, a);

  ARComplexDouble sigma = shift;

  ARluCompStdEig<ARfloat> Prob(numberOfEigenvalues, A, sigma);
  
  if(findEigenvectors)
	Prob.FindEigenvectors();
  else
	Prob.FindEigenvalues();
  
  if(Prob.ConvergedEigenvalues() != (long)numberOfEigenvalues)
	{
	  throw RLException("Only %d eigenvalues converged, %d were requested.\n",Prob.ConvergedEigenvalues(), numberOfEigenvalues);
	}
  EigenInformation * toReturn = new EigenInformation();
  toReturn->Eigenvalues.resize(Prob.ConvergedEigenvalues());
  if(findEigenvectors)
	toReturn->Eigenvectors.resize(Prob.ConvergedEigenvalues(), vector<ComplexDouble>(n, 0.0) );
  for(long i = 0; i<Prob.ConvergedEigenvalues(); ++i)
	{
	  toReturn->Eigenvalues.at(i) = Prob.Eigenvalue(i);
	  if(findEigenvectors)
		for(ulong j = 0; j<n; ++j)
		  {
			toReturn->Eigenvectors.at(i).at(j) = Prob.Eigenvector(i, j);
		  }
	}


  return toReturn;
}
