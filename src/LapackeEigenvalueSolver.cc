#include "LapackeEigenvalueSolver.hh"

RIKARD_COMPLEX_TYPE * LapackeEigenvalueSolver::GetArray(CMatrix * matrix)
{
  ulong N = matrix->Rows() * matrix->Columns();
  RIKARD_COMPLEX_TYPE * toReturn = new RIKARD_COMPLEX_TYPE[N];
  ComplexDouble * a = matrix->GetArray();
  for(ulong i = 0; i<N; ++i)
	{
	  RIKARD_COMPLEX_COPY_TO(a[i], toReturn[i]);
	}
  matrix->DeleteArray();
  return toReturn;
}

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

  RIKARD_COMPLEX_TYPE * a = GetArray(toSolve);

  int reply = LAPACKE_zgeev(matrix_order, jobvl, jobvr,
							n, a, lda, &w[0], vl, ldvl, &vr[0], ldvr);
  
  if( reply )
	{
	  throw RLException("LapackeEigenvalueSolver: LAPACKe terminated with return code %d.", reply);
	}

  EigenInformation * toReturn = new EigenInformation();
  toReturn->Eigenvalues.resize(n);
  for(int i = 0; i<n; ++i)
	{
	  RIKARD_COMPLEX_COPY_FROM(w[i], toReturn->Eigenvalues[i]);
	}
  
  if(ComputeEigenvectors)
	{
	  toReturn->Eigenvectors.resize(n);
	  for(int i = 0; i<n; ++i)
		{
		  toReturn->Eigenvectors[i].resize(n);
		}
	  for(int i = 0; i<n*n; ++i)
		{
		  RIKARD_COMPLEX_COPY_FROM(vr.at(i), toReturn->Eigenvectors[i%n][i/n]);
		}
	}

  delete [] a;

  return toReturn;
}
