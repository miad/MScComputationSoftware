#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Potential.hh"

#include "KPoints.hh"
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


#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);

using namespace std;


ComplexDouble ExpOfInnerProduct(ComplexDouble k1, ComplexDouble k2);

ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2);

int main(int argc, char *argv[]);


CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.


#endif
