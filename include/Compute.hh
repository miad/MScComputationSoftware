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



int main(int argc, char *argv[]);

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.

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
