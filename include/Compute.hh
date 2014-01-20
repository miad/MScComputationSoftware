#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <unistd.h>
#include <algorithm>

#include "Potential.hh"
#include "ParametrizedPotential.hh"
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
#include "HermiteEvaluator.hh"
#include "HarmonicBasisFunction.hh"
#include "CompositeBasisFunction.hh"


using namespace std;


#define TMP_LIST_STR(n, s) list<string> n; n.push_back(s);


int main(int argc, 
		 char *argv[]
		 );


void PerformSolution(ComputeConfig & myConfiguration, 
					 VerbosePrinter & myPrinter, 
					 OutputProcessor & myProcessor
					 );

void VerifyMatrixBasicProperties(ComputeConfig & myConfiguration, 
								 VerbosePrinter & myPrinter, 
								 CMatrix * HamiltonianMatrix
								 );

CommandLineInterpreter * InitInterpreter();///Returns a command line interpreter with defined commands.

void PrintNumberOfNonzeroElements(CMatrix * HamiltonianMatrix, 
								  VerbosePrinter & myPrinter
								  );

void * EvaluateSubMatrixOneParticle(OneParticleWorkerData w
									);

void * EvaluateSubMatrixTwoParticles(TwoParticleWorkerData w
									 );

void * EvaluateSubMatrixOneParticleHarmonic(OneParticleWorkerData w
											);

CMatrix * ConstructOneParticleHamiltonian(const ComputeConfig & myConfiguration, 
										  VerbosePrinter & myPrinter,
										  uint particleID = 0
										  );

CMatrix * ConstructTwoParticleHamiltonian(const ComputeConfig & myConfiguration,
										  VerbosePrinter & myPrinter, 
										  PrecomputedInteractionEvaluator & myPrecomputedInteractionEvaluator
										  
										  );


#endif
