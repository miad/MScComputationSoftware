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
				   ulong _numberOfEigenvalues = 10, 
				   ComplexDouble _shift = 0.0,
				   string _workArea = "", 
				   bool _findEigenvectors = true
				   );

  ~EigenvalueSolver();

  ComplexDouble GetShift() const;

  void SetShift(ComplexDouble value
				);
  
  ulong GetNumberOfEigenvalues() const;

  void SetNumberOfEigenvalues(ulong value
							  );

  EigenSolverType GetSolverType() const;

  void SetSolverType(EigenSolverType value
					 );

  string GetWorkArea() const;
  
  void SetWorkArea(string value
				   );

  bool GetFindEigenvectors() const;

  void SetFindEigenvectors(bool value
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
  ulong numberOfEigenvalues;
  ComplexDouble shift;
  string workArea;
  bool findEigenvectors;
};

#endif
