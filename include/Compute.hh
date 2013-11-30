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



class WorkerData
{
public:
  WorkerData(CMatrix * _HamiltonianMatrix,
			 ParametrizedCurve * _myCurve,
			 Potential * _myPotential,
			 vector<BasisFunction> _myBasisFunctions,
			 unsigned int _numberOfGLPoints,
			 unsigned int _m1, unsigned int _m2,
			 unsigned int _n1, unsigned int _n2)
	:HamiltonianMatrix(_HamiltonianMatrix),
	 myCurve(_myCurve),
	 myPotential(_myPotential),
	 myBasisFunctions(_myBasisFunctions),
	 numberOfGLPoints(_numberOfGLPoints),
	 m1(_m1),m2(_m2),n1(_n1),n2(_n2)
  { 
	for(vector<BasisFunction>::iterator it = myBasisFunctions.begin(); it!=myBasisFunctions.end(); ++it)
	  {
		it->ForceDeepCopy();
	  }
	
}

  CMatrix * HamiltonianMatrix;
  ParametrizedCurve * myCurve; 
  Potential * myPotential;
  vector<BasisFunction> myBasisFunctions;
  unsigned int numberOfGLPoints; 
  unsigned int m1, m2, n1, n2;
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



#endif
