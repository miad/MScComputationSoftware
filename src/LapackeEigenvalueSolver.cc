#include "LapackeEigenvalueSolver.hh"

EigenInformation * LapackeEigenvalueSolver::Solve(CMatrix * toSolve, bool ComputeEigenvectors)
{
  if( ! toSolve->IsSquare() )
	{
	  throw RLException("Can only find eigenvalues of square matrices.");
	}
  
  int matrix_order = LAPACK_ROW_MAJOR;
  lapack_int n = toSolve->Rows();
  ComplexDouble * a = toSolve->GetArray();
  lapack_int lda = n;
  vector<lapack_complex_double> w(n, 0);
  char jobvl = 'N';
  lapack_complex_double * vl = 0;
  lapack_int ldvl = n;
  char jobvr = 'N';
  vector<lapack_complex_double> vr;
  lapack_int ldvr = n;
  
  if(ComputeEigenvectors)
	{
	  jobvr = 'V';
	  vr.resize(n*n, 0);
	}
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
	  for(int i = 0; i<n*n; ++i)
		{
		  if(i/n==0)
			toReturn->Eigenvectors.push_back(vector<ComplexDouble>());
		  toReturn->Eigenvectors[i%n].push_back(vr.at(i));
		}
	}


  return toReturn;
}
