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
#include "SimpsonIntegrator.hpp"
#include <math.h>
#include "EigenvalueSolver.hh"
#include "RLException.hh"
#include "CommandLineInterpreter.hh"
#include "CommandLineArgument.hh"
#include "CommandLineException.hh"
#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"


using namespace std;

Potential * myPotential;


ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2);

int main(int argc, char *argv[]);


CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.

void PrintHelp(); //! Prints help message.


#endif
