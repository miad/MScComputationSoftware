#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Potential.hh"

#include <iostream>
#include "Globals.hpp"
#include "Matrix.hpp"
#include <math.h>
#include "EigenvalueSolver.hh"
#include "RLException.hh"
#include "CommandLineInterpreter.hh"
#include "CommandLineArgument.hh"
#include "CommandLineException.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "ParametrizedCurve.hh"
#include "LegendreRule.hh"
#include <unistd.h>


#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);

#ifdef EPS
#undef EPS
#endif

#define EPS 1E-9

#ifndef PI
#define PI 3.141592653589793238462643
#endif

#define NPPUSH(v) numberOfPointsOnCurve.push_back(v)

#define DBL_EQUAL(a, b) (abs(b-a) < EPS)

using namespace std;


ComplexDouble ExpOfInnerProduct(ComplexDouble k1, 
								ComplexDouble k2
								);

void PrintDataToFile(VerbosePrinter * myPrinter,
					 const string fileName, 
					 const EigenInformation & data, 
					 const vector<ComplexDouble> & kValuesOnCurve,
					 const list<Interval> & potentialIntervals
					 );

bool invertInnerProduct;

int main(int argc, char *argv[]);

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.


#endif
