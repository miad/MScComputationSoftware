#include "Compute.hh"

int main(int argc, char *argv[])
{
  CommandLineInterpreter * myInterpreter = InitInterpreter(); 
  try
	{
	  myInterpreter->Initialize(argc, argv);
	  if(!myInterpreter->ReadFlaggedCommand("help").empty() || argc < 1)
		{
		  myInterpreter->PrintHelp();
		  return 0;
		}
	}
  catch(CommandLineException e)
	{
	  cerr << e.what() << endl;
	  return 1;
	}

  ///Read configuration file and act based on that.

  string configFile = myInterpreter->ReadFlaggedCommandStrict("configFile").front().c_str();

  delete myInterpreter;


  ComputeConfig myConfiguration;
  
  ///If there is no config file, write it and exit. Otherwise, use it.
  if( access(configFile.c_str(), F_OK) == -1 )
	{
	  printf("Config file '%s' not found, creating it and quitting.\n", configFile.c_str());
	  myConfiguration.WriteFile(configFile.c_str());
	  return 0;
	}
  else
	{
	  myConfiguration.ReadFile(configFile.c_str());
	}


  ///Initialize stuff.
  
  VerbosePrinter myPrinter(myConfiguration.GetVerbosityLevel());
  myPrinter.Print(3, "Initializing...\n");
  OutputProcessor myProcessor(&myConfiguration);
  myProcessor.RegisterListener(&myPrinter);


  ///Solve problem and save data.
  PerformSolution(myConfiguration, myPrinter, myProcessor);
  

  ///Launch plotters if we should do so.
  myPrinter.Print(2, "Launching plotters.\n");
  ExternalLauncher myLauncher(&myConfiguration);
  myLauncher.Launch();


  myPrinter.Print(2, "Done, exiting.\n");
  return 0;
}


void PerformSolution(ComputeConfig & myConfiguration, VerbosePrinter & myPrinter, OutputProcessor & myProcessor)
{  
  if(myConfiguration.GetNumberOfParticles() == 1)
	{
	  CMatrix * HamiltonianMatrix = ConstructOneParticleHamiltonian(myConfiguration, myPrinter);
	  VerifyMatrixBasicProperties(myConfiguration, myPrinter, HamiltonianMatrix);

	  PrintNumberOfNonzeroElements(HamiltonianMatrix, myPrinter);
	  myProcessor.SaveMatrix(HamiltonianMatrix);
	  
	  myPrinter.Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");
	  EigenInformation * myInfo = myConfiguration.GetSolver()->Solve(HamiltonianMatrix);
	  
	  myPrinter.Print(2, "Saving results to files.\n");
	  
	  myProcessor.SetEigenInformation(myInfo);
	  myProcessor.WritePostOutput();
	  delete HamiltonianMatrix;
	  delete myInfo;
	}
  else if(myConfiguration.GetNumberOfParticles() == 2)
	{
	  PrecomputedInteractionEvaluator *  myPrecomputedInteractionEvaluator = NULL;
	  vector<CompositeBasisFunction * > myBasisFunctions;
	  vector<EigenInformation *> OneParticleEigenData;
	  int r[2]; r[0] = -1; r[1] = -1;
	  
	  if( myConfiguration.GetInteractionProperties()->GetCacheFile().empty() ||
		  access(myConfiguration.GetInteractionProperties()->GetCacheFile().c_str(), F_OK) == -1)
		{
		  for(uint i = 0; i<2; ++i)
			{
			  myPrinter.Print(1, "Constructing Hamiltonian for particle %d.\n", i);
			  CMatrix * OneBodyHamiltonian = ConstructOneParticleHamiltonian(myConfiguration, myPrinter, i);
			  VerifyMatrixBasicProperties(myConfiguration, myPrinter, OneBodyHamiltonian);

			  myPrinter.Print(1, "Finding eigenvalues for particle %d.\n", i);
			  EigenInformation * myEigenInfo = myConfiguration.GetSolver()->LapackeSolve(OneBodyHamiltonian);
			  if(myConfiguration.GetHarmonicOverride())
				{
				  myBasisFunctions.push_back(new CompositeBasisFunction(
																		myConfiguration.GetHarmonicBasisFunction(),
																		myEigenInfo,
																		i
																		)
											 );
				}
			  else
				{
				  myBasisFunctions.push_back(new CompositeBasisFunction(
																		myConfiguration.GetBasisFunctions(),
																		myEigenInfo,
																		myConfiguration.GetKCurve()
																		)
											 );
				}
			  myProcessor.SetEigenInformation(myEigenInfo);
			  myProcessor.WriteInterestingKPointsVerbosely(true, i);
			  vector<uint> idx = myProcessor.FindInterestingKPointIndex(true);
			  if(idx.empty())
				{
				  myPrinter.Print(2, "Warning: no interesting K point for further processing.\n");
				  r[0] = -10;
				}
			  else
				{
				  if(idx.size() > 1)
					myPrinter.Print(2, "Warning: too many interesting K points, selecting first one.\n");
				  r[i] = idx.front();
				}
			  
			  OneParticleEigenData.push_back(myEigenInfo);  ///Must be retained in memory and deleted later.
			  delete OneBodyHamiltonian;
			  
			}
		


	  
		  myPrinter.Print(2, "Precomputing interaction terms:\n");
		  myPrecomputedInteractionEvaluator = 
			new PrecomputedInteractionEvaluator(&myBasisFunctions, myConfiguration.GetInteractionProperties(), myConfiguration.GetHarmonicBasisFunction(), myConfiguration.GetSpecificUnits(), &myPrinter, myConfiguration.GetNumberOfThreads());
		  
		  myPrinter.Print(2, "Done precomputing.\n");
		}
	  else
		{
		  myPrinter.Print(2, "Loading interaction properties directly from save file.\n");
		  myPrecomputedInteractionEvaluator = new PrecomputedInteractionEvaluator(myConfiguration.GetInteractionProperties(), &myPrinter);
		  myPrinter.Print(2, "Done loading interaction properties.\n");
		}

	  CMatrix * TwoBodyHamiltonian = ConstructTwoParticleHamiltonian(myConfiguration, myPrinter, myPrecomputedInteractionEvaluator);
	  VerifyMatrixBasicProperties(myConfiguration, myPrinter, TwoBodyHamiltonian);

	  ///Free some memory before invoking eigenvalue solver.
	  delete myPrecomputedInteractionEvaluator; myPrecomputedInteractionEvaluator = NULL;


	  myPrinter.Print(2, "Optionally saving Lanczos object...");
	  ///If we should save and exit, do so.
	  bool doBreak = false;
	  doBreak = LanczosSaver::Save(myConfiguration.GetOutputFilenames()->Get("LanczosObject").c_str(),
								  TwoBodyHamiltonian,
								   r
								  );
	  myPrinter.Print(2, "done.\n");
	  doBreak = myProcessor.SaveMatrix(TwoBodyHamiltonian) || doBreak; 

	  if(doBreak)
		return;



	  myPrinter.Print(1, "Finding eigenvalues for the two-body system.\n");

	  ///NOTE: I have disabled the throw() on non-orthogonal eigenvectors because the eigenvectors are NOT orthogonal any more...
	  EigenInformation * TwoBodyEigenInfo = myConfiguration.GetSolver()->Solve(TwoBodyHamiltonian, false);

	  myConfiguration.GetSolver()->RescaleEigenvectors(TwoBodyEigenInfo);

	  myProcessor.SetEigenInformation(TwoBodyEigenInfo);

	  myProcessor.SetCompositeBasisFunctions(&myBasisFunctions);




	  myProcessor.WritePostOutput();
	  

	  /// finally...
	  delete TwoBodyHamiltonian;
	  delete TwoBodyEigenInfo;
	  for(vector<EigenInformation * >::iterator it = OneParticleEigenData.begin(); it!=OneParticleEigenData.end(); ++it)
		delete *it;

	}
  else
	{
	  throw RLException("Invalid number of particles: %d", myConfiguration.GetNumberOfParticles());
	}
}


void VerifyMatrixBasicProperties(ComputeConfig & myConfiguration, VerbosePrinter & myPrinter, CMatrix * HamiltonianMatrix)
{
  ///Validate some basic properties of the Hamilton matrix, if they are expected.
  if(myConfiguration.GetExpectedMatrixType() == SymmetricMatrix)
	{
	  myPrinter.Print(1, "Validating symmetricity of matrix.\n");
	  if ( ! HamiltonianMatrix->IsSymmetric(true) )
		{
		  throw RLException("The matrix was found to be non-symmetric.");
		}
	}
  if(myConfiguration.GetExpectedMatrixType() == HermitianMatrix)
	{
	  myPrinter.Print(1, "Validating hermiticity of matrix.\n");
	  if( ! HamiltonianMatrix->IsHermitian(true))
		{
		  throw RLException("The matrix was found to be non-hermitian.");
		}
	}

  if(myConfiguration.GetHarmonicOverride())
	{
	  myPrinter.Print(1, "Harmonic: Validating reality of matrix.\n");
	  for(uint i = 0; i<HamiltonianMatrix->Rows(); ++i)
		{
		  for(uint j = 0; j<HamiltonianMatrix->Columns(); ++j)
			{
			  if(imag(HamiltonianMatrix->Element(i, j)) > 1E-9)
				throw RLException("Element (%d, %d) had an imaginary part %f, but a real matrix was expected.", i, j, imag(HamiltonianMatrix->Element(i, j)));
			}
		}
	}


}



void PrintNumberOfNonzeroElements(CMatrix * HamiltonianMatrix, VerbosePrinter & myPrinter)
{
  uint nonzero = 0;
  for(uint i = 0; i<HamiltonianMatrix->Rows(); ++i)
	{
	  for(uint j = 0; j<HamiltonianMatrix->Columns(); ++j)
		{
		  nonzero += !DBL_EQUAL(HamiltonianMatrix->Element(i, j),0.0);
		}
	}
  myPrinter.Print(3, "Number of nonzero Hamiltonian elements: %d\n", nonzero);
}




CMatrix * ConstructOneParticleHamiltonian(const ComputeConfig & myConfiguration, VerbosePrinter & myPrinter, uint particleID)
{

  vector<BasisFunction> myBasisFunctions = myConfiguration.GetBasisFunctions();
  if(myConfiguration.GetHarmonicOverride())
	myPrinter.Print(1, "Harmonic override enabled, some functions will alter behavior.\n");
  else
	myPrinter.Print(5, "Using %d basis functions:", myBasisFunctions.size());

  for(vector<BasisFunction>::const_iterator it = myBasisFunctions.begin(); it!=myBasisFunctions.end(); ++it)
	{
	  myPrinter.Print(5, "%s %s", ((it==myBasisFunctions.begin())?"":","),it->GetName());
	}
  myPrinter.Print(5, ".\n");

  if(particleID >= myConfiguration.GetNumberOfParticles())
	{
	  throw RLException("ConstructOneParticleHamiltonian called with invalid particle ID: %d", particleID);
	}

  myPrinter.Print(3, "Creating 1-particle Hamiltonian for particle : %d\n", particleID);
  
  uint numberOfGLPoints = myConfiguration.GetKCurve()->GetTotalNumberOfGLPoints();
  uint numberOfBasisFunctions =  (uint)myBasisFunctions.size();
  uint MatrixSize = numberOfGLPoints * numberOfBasisFunctions;

  if(myConfiguration.GetHarmonicOverride())
	MatrixSize = myConfiguration.GetHarmonicNmax() ;

  myPrinter.Print(5, "Total number of GL points in k-space: %d, Matrix size: %d.\n", numberOfGLPoints, MatrixSize);


  CMatrix * HamiltonianMatrix = new CMatrix(MatrixSize, MatrixSize);
  HamiltonianMatrix->InitializeAll(0.);


  ///The main stuff: construct the Hamiltoninan matrix.

  myPrinter.Print(1,"Constructing Hamiltonian matrix.\n");
  

  ParametrizedCurve * myCurve = myConfiguration.GetKCurve();
  Potential * myPotential = myConfiguration.GetPotential(particleID);

  myPrinter.Print(3, "Potential: minX %13.10e maxX %13.10e\n", myPotential->GetMinX(), myPotential->GetMaxX());

  MultiTasker<OneParticleWorkerData, void*> * myMultiTasker = NULL;
  ///Generate the Hamiltonian in parallell!
  if(myConfiguration.GetHarmonicOverride())
	{
	  HermiteEvaluator::Init(myConfiguration.GetHarmonicNmax() + 10);
	  myMultiTasker = new MultiTasker<OneParticleWorkerData, void*>(EvaluateSubMatrixOneParticleHarmonic, myConfiguration.GetNumberOfThreads());
	}
  else
	{
	  myMultiTasker = new MultiTasker<OneParticleWorkerData, void*>(EvaluateSubMatrixOneParticle, myConfiguration.GetNumberOfThreads());
	}

  if(myMultiTasker == NULL)
	{
	  throw RLException("This should never happen.");
	}
  
  myMultiTasker->RegisterListener(&myPrinter);
  
  for(uint i = 0; i<MatrixSize; ++i)
	{
	  myMultiTasker->AddInput(OneParticleWorkerData(HamiltonianMatrix,
													myCurve,
													myPotential,
													myBasisFunctions,
													numberOfGLPoints,
													myConfiguration.GetSpecificUnits()->GetHbarTimesLambda(),
													myConfiguration.GetSpecificUnits()->GetMassOverLambda2(),
													myConfiguration.GetHarmonicBasisFunction(),
													particleID,
													i,
													i+1,
													0,
													MatrixSize
													)
							  );
	  
	}
  myMultiTasker->LaunchThreads();
  myMultiTasker->PauseUntilOutputIsGenerated();
  myMultiTasker->DestroyThreads();
  delete myMultiTasker; 
  myMultiTasker = NULL;

  return HamiltonianMatrix;
}




CMatrix * ConstructTwoParticleHamiltonian(const ComputeConfig & myConfiguration, VerbosePrinter & myPrinter, PrecomputedInteractionEvaluator * myPrecomputedInteractionEvaluator)
{

  myPrinter.Print(3, "Creating 2-particle Hamiltonian\n");
  
  uint MatrixSize = myPrecomputedInteractionEvaluator->GetBasisElements(0) * myPrecomputedInteractionEvaluator->GetBasisElements(1);

  CMatrix * HamiltonianMatrix = new CMatrix(MatrixSize, MatrixSize);
  HamiltonianMatrix->InitializeAll(0.);



  MultiTasker<TwoParticleWorkerData, void*> * myMultiTasker = NULL;

  ///Generate the Hamiltonian in parallell!
  myMultiTasker = new MultiTasker<TwoParticleWorkerData, void*>(EvaluateSubMatrixTwoParticles, myConfiguration.GetNumberOfThreads());

  myMultiTasker->RegisterListener(&myPrinter);
  

  for(uint i = 0; i<MatrixSize; ++i)
	{
	  myMultiTasker->AddInput(TwoParticleWorkerData(HamiltonianMatrix,
													myPrecomputedInteractionEvaluator,
													i,
													i+1,
													0,
													MatrixSize
													)
							  );
	  
	}
  myMultiTasker->LaunchThreads();
  myMultiTasker->PauseUntilOutputIsGenerated();
  myMultiTasker->DestroyThreads();
  delete myMultiTasker; 
  myMultiTasker = NULL;


  myPrinter.Print(3, "Done, returning Hamiltonian.\n");
  return HamiltonianMatrix;
}





CommandLineInterpreter * InitInterpreter()
{
  CommandLineInterpreter * myInterpreter = new CommandLineInterpreter();

  TMP_LIST_STR(confDefault,"config.conf");
  TMP_LIST_STR(confDescription,"Configuration file");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("configFile", 1, false, "Location of program config file. If the file does not exist, it is created with a default configuration.", confDescription, confDefault));

  myInterpreter->AddCommandLineArgument(CommandLineArgument("help",0,false, "Displays a help message and quits."));

  myInterpreter->SetDescription("Compute eigenvalues for particle states in a potential.");
  return myInterpreter;
}

void * EvaluateSubMatrixOneParticle(OneParticleWorkerData w)
{
  CMatrix * HamiltonianMatrix = w.HamiltonianMatrix;
  ParametrizedCurve * myCurve = w.myCurve;
  Potential * myPotential = w.myPotential;
  vector<BasisFunction> * myBasisFunctions = &w.myBasisFunctions;
  uint numberOfGLPoints = w.numberOfGLPoints;
  double hbarTimesLambda = w.hbarTimesLambda;
  double massOverLambda2 = w.massOverLambda2;
  uint m1 = w.m1, m2 = w.m2, n1 = w.n1, n2 = w.n2;

  for(uint i = m1; i<m2; ++i)
	{
	  uint curvePointerA = i % numberOfGLPoints;
	  uint basisPointerA = i / numberOfGLPoints;
	  uint curveSegmentA = myCurve->SegmentIndexFromGLNumber(curvePointerA);
	  
	  ComplexDouble kA = myCurve->GetRuleValue(curveSegmentA, curvePointerA);
	  ComplexDouble wA = myCurve->GetRuleWeight(curveSegmentA, curvePointerA);
	  
	  for(uint j = n1; j<n2; ++j)
		{
		  uint curvePointerB = j % numberOfGLPoints;
		  uint basisPointerB = j / numberOfGLPoints;
		  uint curveSegmentB = myCurve->SegmentIndexFromGLNumber(curvePointerB);
		  
		  ComplexDouble kB = myCurve->GetRuleValue(curveSegmentB, curvePointerB);
		  ComplexDouble wB = myCurve->GetRuleWeight(curveSegmentB, curvePointerB);
		  
		  HamiltonianMatrix->Element(i, j) += sqrt(wA*wB)*
			myPotential->BasisIntegrate((*myBasisFunctions)[basisPointerA], 
										(*myBasisFunctions)[basisPointerB], 
										kA, kB) ;
		}
	HamiltonianMatrix->Element(i,i) += pow(hbarTimesLambda,2)/(2.*massOverLambda2) * pow(kA, 2);
   }
   return NULL;
}


void * EvaluateSubMatrixTwoParticles(TwoParticleWorkerData w)
{
  CMatrix * HamiltonianMatrix = w.HamiltonianMatrix;
  const PrecomputedInteractionEvaluator * myPrecomputedInteractionEvaluator = w.myPrecomputedInteractionEvaluator;
  uint m1 = w.m1, m2 = w.m2, n1 = w.n1, n2 = w.n2;

  if(myPrecomputedInteractionEvaluator == NULL)
	{
	  throw RLException("myPrecomputedInteractionEvaluator was NULL.");
	}

  uint N1 = myPrecomputedInteractionEvaluator->GetBasisElements(0);
  uint N2 = myPrecomputedInteractionEvaluator->GetBasisElements(1);


  for(uint i = m1; i<m2; ++i)
	{
	  uint a = i / N1;
	  uint b = i % N2;

	  for(uint j = n1; j<n2; ++j)
		{
		  uint c = j / N1;
		  uint d = j % N2;
		  
		  if(a==c && b==d)
			{
			  HamiltonianMatrix->Element(i, j) += myPrecomputedInteractionEvaluator->GetEnergy(0, a);
			  HamiltonianMatrix->Element(i, j) += myPrecomputedInteractionEvaluator->GetEnergy(1, b);
			}
		  
		  ///Here we have to be very careful not to cause any double-addings. 
		  HamiltonianMatrix->Element(i, j) += myPrecomputedInteractionEvaluator->GetElement(a, b, c, d);
			  
		}
	}
  return NULL;
}


void * EvaluateSubMatrixOneParticleHarmonic(OneParticleWorkerData w)
{
  CMatrix * HamiltonianMatrix = w.HamiltonianMatrix;

  uint m1 = w.m1, m2 = w.m2, n1 = w.n1, n2 = w.n2;
  for(uint i = m1; i<m2; ++i)
	{
	  for(uint j = n1; j<n2; ++j)
		{
		  ///Kinetic part.
		  HamiltonianMatrix->Element(i, j) += w.myHarmonicBasisFunction->KineticTerm(i, j, w.particleID);

		  ///Potential part.
		  HamiltonianMatrix->Element(i, j) += w.myHarmonicBasisFunction->Integrate(i, j, w.particleID);
		}
   }
  return NULL;
}
