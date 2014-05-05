#include "EigenvalueSolver.hh"

EigenvalueSolver::EigenvalueSolver(EigenSolverType _mySolverType, ulong _numberOfEigenvalues, ComplexDouble _shift, string _workArea, bool _findEigenvectors)
  :mySolverType(_mySolverType), numberOfEigenvalues(_numberOfEigenvalues), shift(_shift), workArea(_workArea)
{

}


EigenvalueSolver::~EigenvalueSolver()
{

}

ComplexDouble EigenvalueSolver::GetShift() const
{
  return shift;
}

void EigenvalueSolver::SetShift(ComplexDouble value)
{
  shift = value;
}


ulong EigenvalueSolver::GetNumberOfEigenvalues() const
{
  return numberOfEigenvalues;
}


void EigenvalueSolver::SetNumberOfEigenvalues(ulong value)
{
  numberOfEigenvalues = value;
}

EigenSolverType EigenvalueSolver::GetSolverType() const
{
  return mySolverType;
}

void EigenvalueSolver::SetSolverType(EigenSolverType value)
{
  mySolverType = value;
}

string EigenvalueSolver::GetWorkArea() const
{
  return workArea;
}

void EigenvalueSolver::SetWorkArea(string value)
{
  workArea = value;
}

EigenInformation * EigenvalueSolver::LapackeSolve(CMatrix * toSolve, bool assureEigenOrthonormality) const
{
  EigenInformation * toReturn = LapackeEigenvalueSolver::Solve(toSolve);

  if(assureEigenOrthonormality)
	{
	  AssureEigenOrthonormality(toReturn);
	}
  return toReturn;

}

EigenInformation * EigenvalueSolver::Solve(CMatrix * toSolve, bool assureEigenOrthonormality) const
{
  EigenInformation * toReturn = NULL;
  switch(mySolverType)
	{
	case LapackSolver:
	  {
		toReturn = LapackeEigenvalueSolver::Solve(toSolve, findEigenvectors);
	  }
	  break;
	case ArpackSolver:
	  {
		toReturn = ArpackEigenvalueSolver::Solve(toSolve, numberOfEigenvalues, shift, findEigenvectors);
	  }
	  break;
	default:
	  {
		throw RLException("Unsupported eigensolver type.");
	  }
	}


  if(assureEigenOrthonormality)
	{
	  AssureEigenOrthonormality(toReturn);
	}
  return toReturn;
}


bool EigenvalueSolver::GetFindEigenvectors() const
{
  return findEigenvectors;
}

void EigenvalueSolver::SetFindEigenvectors(bool value)
{
  findEigenvectors = value;
}



void EigenvalueSolver::RescaleEigenvectors(EigenInformation * eigenData)
{
  ///First rescale all eigenvectors to norm 1.
  for(ulong i = 0; i<eigenData->EigenPairs.size(); ++i)
	{
	  ComplexDouble sum = 0.0;
	  for(ulong a = 0; a<eigenData->EigenPairs[i].Eigenvector.size(); ++a)
		{
		  sum += pow(eigenData->EigenPairs[i].Eigenvector[a], 2);
		}
	  ComplexDouble sqRoot = sqrt(sum);
	  for(ulong a = 0; a<eigenData->EigenPairs[i].Eigenvector.size(); ++a)
		{
		  eigenData->EigenPairs[i].Eigenvector[a] /= sqRoot;
		}	  
	}
}


double EigenvalueSolver::AssureEigenOrthonormality(EigenInformation * eigenData)
{
  RescaleEigenvectors(eigenData);

  ///Used for exception info on failure.
  double maxDeviation = 0;
  pair<ulong, ulong> maxP;

  ///Now verify orthonormality and throw if fail.
  for(ulong i = 0; i<eigenData->EigenPairs.size(); ++i)
	{
	  for(ulong j = i; j<eigenData->EigenPairs.size(); ++j)
		{
		  ComplexDouble sum = 0.0;
		  for(ulong a = 0; a<eigenData->EigenPairs[i].Eigenvector.size(); ++a)
			{
			  sum += eigenData->EigenPairs[i].Eigenvector[a] * eigenData->EigenPairs[j].Eigenvector[a];
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
	  throw RLException("Too large eigenvector deviation: %+13.10e between eigenvectors %d and %d with eigenvalues %+13.10f%+13.10fi and %+13.10f%+13.10fi respectively. Investigate this.", maxDeviation, maxP.first, maxP.second, real(eigenData->EigenPairs[maxP.first].Eigenvalue), imag(eigenData->EigenPairs[maxP.first].Eigenvalue), real(eigenData->EigenPairs[maxP.second].Eigenvalue), imag(eigenData->EigenPairs[maxP.second].Eigenvalue));
	}
  return maxDeviation;
}
