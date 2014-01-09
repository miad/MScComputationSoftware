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
#include "LapackeEigenvalueSolver.hh"
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
#include "WorkerData.hh"
#include "OutputProcessor.hh"
#include "SpecificUnits.hh"
#include "ExternalLauncher.hh"


using namespace std;


#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);

#define DELTA(a, b) ((double)(a==b))


int main(int argc, char *argv[]);

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.

void * EvaluateSubMatrixOneParticle(WorkerData w);

void * EvaluateSubMatrixTwoParticles(WorkerData w);

void * EvaluateSubMatrixOneParticleHarmonic(WorkerData w);


CMatrix * ConstructHamiltonian(const ComputeConfig & myConfiguration, 
							   VerbosePrinter & myPrinter);


#endif
