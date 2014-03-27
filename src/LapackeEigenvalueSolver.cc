#include "LapackeEigenvalueSolver.hh"

EigenInformation * LapackeEigenvalueSolver::Solve(CMatrix * toSolve, bool ComputeEigenvectors)
{
  if( ! toSolve->IsSquare() )
	{
	  throw RLException("Can only find eigenvalues of square matrices.");
	}
  
  int matrix_order = LAPACK_ROW_MAJOR;
  lapack_int n = toSolve->Rows();
  lapack_int lda = n;
  vector<RIKARD_COMPLEX_TYPE> w(n);
  char jobvl = 'N';
  RIKARD_COMPLEX_TYPE * vl = 0;
  lapack_int ldvl = n;
  char jobvr = 'N';
  vector<RIKARD_COMPLEX_TYPE> vr;
  lapack_int ldvr = n;
  
  if(ComputeEigenvectors)
	{
	  jobvr = 'V';
	  vr.resize(n*n);
	}

  RIKARD_COMPLEX_TYPE * a = toSolve->GetArray();

  int reply = LAPACKE_zgeev(matrix_order, jobvl, jobvr,
							n, a, lda, &w[0], vl, ldvl, &vr[0], ldvr);
  
  if( reply )
	{
	  throw RLException("LapackeEigenvalueSolver: LAPACKe terminated with return code %d.", reply);
	}

  EigenInformation * toReturn = new EigenInformation();
  toReturn->Eigenvalues = w;

  if(ComputeEigenvectors)
	{
	  toReturn->Eigenvectors.resize(n);
	  for(uint i = 0; i<n; ++i)
		{
		  toReturn->Eigenvectors[i].resize(n);
		}
	  for(uint i = 0; i<n*n; ++i)
		{
		  toReturn->Eigenvectors[i%n][i/n] = vr.at(i);
		}
	}

  return toReturn;
}
