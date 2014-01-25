#include "LapackeEigenvalueSolver.hh"

EigenInformation * LapackeEigenvalueSolver::Solve(CMatrix * toSolve, bool assureEigenOrthonormality)
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
  char jobvr = 'V';
  vector<lapack_complex_double> vr(n*n,0);
  lapack_int ldvr = n;
  
  int reply = LAPACKE_zgeev(matrix_order, jobvl, jobvr,
							n, a, lda, &w[0], vl, ldvl, &vr[0], ldvr);
  
  if( reply )
	{
	  throw RLException("LapackeEigenvalueSolver: LAPACKe terminated with return code %d.", reply);
	}

  EigenInformation * toReturn = new EigenInformation();
  toReturn->Eigenvalues = w;
  for(int i = 0; i<n*n; ++i)
	{
	  if(i/n==0)
		toReturn->Eigenvectors.push_back(vector<ComplexDouble>());
	  toReturn->Eigenvectors[i%n].push_back(vr.at(i));
	}

  if(assureEigenOrthonormality)
	{
	  AssureEigenOrthonormality(toReturn);
	}


  return toReturn;
}

void LapackeEigenvalueSolver::RescaleEigenvectors(EigenInformation * eigenData)
{
  ///First rescale all eigenvectors to norm 1.
  for(uint i = 0; i<eigenData->Eigenvectors.size(); ++i)
	{
	  ComplexDouble sum = 0.0;
	  for(uint a = 0; a<eigenData->Eigenvectors[i].size(); ++a)
		{
		  sum += pow(eigenData->Eigenvectors[i][a], 2);
		}
	  ComplexDouble sqRoot = sqrt(sum);
	  for(uint a = 0; a<eigenData->Eigenvectors[i].size(); ++a)
		{
		  eigenData->Eigenvectors[i][a] /= sqRoot;
		}	  
	}
}


double LapackeEigenvalueSolver::AssureEigenOrthonormality(EigenInformation * eigenData)
{
  RescaleEigenvectors(eigenData);


  double maxDeviation = 0;
  pair<uint, uint> maxP;
  ///Now verify orthonormality and throw if fail.
  for(uint i = 0; i<eigenData->Eigenvectors.size(); ++i)
	{
	  for(uint j = i; j<eigenData->Eigenvectors.size(); ++j)
		{
		  ComplexDouble sum = 0.0;
		  for(uint a = 0; a<eigenData->Eigenvectors[i].size(); ++a)
			{
			  sum += eigenData->Eigenvectors[i][a] * eigenData->Eigenvectors[j][a]; 
			}
		  double d = abs(sum-(double)(i==j));
		  if(d > maxDeviation)
			{
			  maxDeviation = d;
			  maxP = make_pair(i, j);
			}
		}
	}
  if(maxDeviation > 1E-6) ///If it's worse than this, it's probably so bad that we should throw() on it.
	{
	  throw RLException("Too large eigenvector deviation: %+13.10e between eigenvectors %d and %d with eigenvalues %+13.10f%+13.10fi and %+13.10f%+13.10fi respectively. Investigate this.", maxDeviation, maxP.first, maxP.second, real(eigenData->Eigenvalues[maxP.first]), imag(eigenData->Eigenvalues[maxP.first]), real(eigenData->Eigenvalues[maxP.second]), imag(eigenData->Eigenvalues[maxP.second]));
	}
  return maxDeviation;
}
