#include "EigenvalueSolver.hh"

EigenvalueSolver::EigenvalueSolver(EigenSolverType _mySolverType, uint _numberOfEigenvalues, ComplexDouble _shift, string _workArea, bool _findEigenvectors)
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


uint EigenvalueSolver::GetNumberOfEigenvalues() const
{
  return numberOfEigenvalues;
}


void EigenvalueSolver::SetNumberOfEigenvalues(uint value)
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


double EigenvalueSolver::AssureEigenOrthonormality(EigenInformation * eigenData)
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
