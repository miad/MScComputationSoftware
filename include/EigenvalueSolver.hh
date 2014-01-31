#ifndef EigenvalueSolver_hh
#define EigenvalueSolver_hh 1

#include "LapackeEigenvalueSolver.hh"
#include "ArpackEigenvalueSolver.hh"
#include "RLMacros.hpp"
#include "Matrix.hpp"
#include "EigenInformation.hh"


enum EigenSolverType
  {
	LapackSolver = 0,
	ArpackSolver = 1
  };

class EigenvalueSolver
{
public:
  EigenvalueSolver(EigenSolverType _mySolverType = LapackSolver, 
				   uint _numberOfEigenvalues = 10, 
				   ComplexDouble _shift = 0.0,
				   string _workArea = ""
				   );

  ~EigenvalueSolver();

  ComplexDouble GetShift() const;

  void SetShift(ComplexDouble value
				);
  
  uint GetNumberOfEigenvalues() const;

  void SetNumberOfEigenvalues(uint value
							  );

  EigenSolverType GetSolverType() const;

  void SetSolverType(EigenSolverType value
					 );

  string GetWorkArea() const;
  
  void SetWorkArea(string value
				   );

  EigenInformation * Solve(CMatrix * toSolve, 
						   bool assureEigenOrthonormality = true
						   ) const;

  EigenInformation * LapackeSolve(CMatrix * toSolve,
						   bool assureEigenOrthoNormality = true
						   ) const;

  static void RescaleEigenvectors(EigenInformation * eigenData
								  );

  static double AssureEigenOrthonormality(EigenInformation * eigenData
										);

private:
  EigenSolverType mySolverType;
  uint numberOfEigenvalues;
  ComplexDouble shift;
  string workArea;
};

#endif
