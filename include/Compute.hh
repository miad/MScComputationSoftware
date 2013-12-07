#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <unistd.h>

#include "Potential.hh"
#include "Globals.hpp"
#include "Matrix.hpp"
#include "EigenvalueSolver.hh"
#include "RLException.hh"
#include "MultiTasker.hpp"
#include "CommandLineInterpreter.hh"
#include "CommandLineArgument.hh"
#include "CommandLineException.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "ParametrizedCurve.hh"
#include "LegendreRule.hh"
#include "ComputeConfig.hh"


using namespace std;


#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);


/**
Class containing the data used by each worker thread.
Essential for the use in RLLib.
*/
class WorkerData
{
public:
  WorkerData(CMatrix * _HamiltonianMatrix, ///Pointer to the Hamilton matrix in use.
			 ParametrizedCurve * _myCurve, ///Pointer to the ParametrizedCurve in use.
			 Potential * _myPotential, ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
			 vector<BasisFunction> _myBasisFunctions, ///Basis functions in use.
			 unsigned int _numberOfGLPoints, ///Number of GL points in use.
			 unsigned int _m1, /// Specifies the submatrix to use.
			 unsigned int _m2, /// Specifies the submatrix to use.
			 unsigned int _n1, /// Specifies the submatrix to use.
			 unsigned int _n2 /// Specifies the submatrix to use.
			 ) ///Constructor, basically initialize all values.
	:HamiltonianMatrix(_HamiltonianMatrix),
	 myCurve(_myCurve),
	 myPotential(_myPotential),
	 myBasisFunctions(_myBasisFunctions),
	 numberOfGLPoints(_numberOfGLPoints),
	 m1(_m1),m2(_m2),n1(_n1),n2(_n2)
  { }
  
  CMatrix * HamiltonianMatrix; ///Pointer to the Hamilton matrix in use.
  ParametrizedCurve * myCurve; ///Pointer to the ParametrizedCurve in use.
  Potential * myPotential; ///Pointer to the potential in use. NOTE: may not be threadsafe, should be checked (due to underlying function evaluator class).
  vector<BasisFunction> myBasisFunctions;
  unsigned int numberOfGLPoints; ///Number of GL points in use.
  unsigned int m1; /// Specifies the submatrix to use.
  unsigned int m2; /// Specifies the submatrix to use.
  unsigned int n1; /// Specifies the submatrix to use.
  unsigned int n2; /// Specifies the submatrix to use.
};




int main(int argc, char *argv[]);

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.


void * EvaluateSubMatrix(WorkerData w);

///TODO: Fix this function. 
void PrintDataToFile(VerbosePrinter * myPrinter,
					 const string fileName, 
					 const EigenInformation & data, 
					 const vector<ComplexDouble> & kValuesOnCurve,
					 const list<Interval> & potentialIntervals
					 );

void PrintPotentialToFile(const char * fileName, 
						  const Potential * potential
						  );

void PrintPotentialPrecisionToFile(const char * fileName, 
						  const Potential * potential
						  );

void PrintParametrizedCurveToFile(const char * fileName, 
								  const ParametrizedCurve * toPrint
								  );

void PrintKCurveToFile(const char * fileName, 
					   const vector<ComplexDouble> & toPrint
					   );

void PrintKFoundToFile(const char * fileName, 
					   const EigenInformation * toPrint
					   );


void PrintInterestingKPointsVerbosely(VerbosePrinter * printer,
									  const EigenInformation * toPrint, 
									  const ParametrizedCurve * filter
									  );

void PrintInterestingKPointsToFile(const char * fileName, 
								   const EigenInformation * toPrint, 
								   const ParametrizedCurve * filter,
								   unsigned int numberOfBasisVectors
								   );

vector<ComplexDouble> FindInterestingKPoints(const EigenInformation * found, 
											 const ParametrizedCurve * filter
											 ); ///Find the interesting k-points (bound states or resonant states) by applying a filter on all eigenvalues.

void PrintInterestingWavefunctionsToFile(VerbosePrinter * printer,
										 const char * fileName,
										 const EigenInformation * toPrint,
										 const ParametrizedCurve * filter,
										 const vector<BasisFunction> * myBasisFunctions,
										 double minX,
										 double maxX,
										 double deltaX
										 );

vector<double> GetBasisRatio(unsigned int numberOfBasisVectors, 
							 const vector<ComplexDouble> toSum, 
							 double & totalSum
							 );


#endif
