#include "Compute.hh"

int main(int argc, char *argv[])
{
  
  ///Initialization: read command line and act based on that.
  
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




  vector<BasisFunction> myBasisFunctions = myConfiguration.GetBasisFunctions();
  myPrinter.Print(5, "Using %d basis functions:", myBasisFunctions.size());
  for(vector<BasisFunction>::const_iterator it = myBasisFunctions.begin(); it!=myBasisFunctions.end(); ++it)
	{
	  myPrinter.Print(5, "%s %s", ((it==myBasisFunctions.begin())?"":","),it->GetName());
	}
  myPrinter.Print(5, ".\n");


  CMatrix * HamiltonianMatrix = ConstructHamiltonian(myConfiguration, myPrinter);

  uint nonzero = 0;
  for(uint i = 0; i<HamiltonianMatrix->Rows(); ++i)
	{
	  for(uint j = 0; j<HamiltonianMatrix->Columns(); ++j)
		{
		  nonzero += !DBL_EQUAL(HamiltonianMatrix->Element(i, j),0.0);
		}
	}
  myPrinter.Print(3, "Number of nonzero Hamiltonian elements: %d\n", nonzero);

  OutputProcessor myProcessor(&myConfiguration);
  myProcessor.RegisterListener(&myPrinter);
  myProcessor.SaveMatrix(HamiltonianMatrix);
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

  myPrinter.Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");
  EigenInformation myInfo = LapackeEigenvalueSolver::Solve(HamiltonianMatrix);
  
  
  myPrinter.Print(2, "Saving results to files.\n");

  myProcessor.SetEigenInformation(&myInfo);
  myProcessor.WritePostOutput();

  delete HamiltonianMatrix;

  myPrinter.Print(2, "Launching plotters.\n");
  
  ExternalLauncher myLauncher(&myConfiguration);
  myLauncher.Launch();

  myPrinter.Print(2, "Done, exiting.\n");
  return 0;
}


CMatrix * ConstructHamiltonian(const ComputeConfig & myConfiguration, VerbosePrinter & myPrinter)
{
  uint numberOfParticles = myConfiguration.GetNumberOfParticles();
  if(numberOfParticles < 1 || numberOfParticles > 2)
	throw RLException("Invalid number of particles.");
  myPrinter.Print(3, "Number of particles: %d\n", numberOfParticles);
  
  vector<BasisFunction> myBasisFunctions = myConfiguration.GetBasisFunctions();

  uint numberOfGLPoints = myConfiguration.GetKCurve()->GetTotalNumberOfGLPoints();
  uint numberOfBasisFunctions =  myBasisFunctions.size();
  uint MatrixSize = numberOfGLPoints * numberOfBasisFunctions;

  if(numberOfParticles == 2)
	MatrixSize *= MatrixSize;

  myPrinter.Print(5, "Total number of GL points in k-space: %d, Matrix size: %d.\n", numberOfGLPoints, MatrixSize);


  CMatrix * HamiltonianMatrix = new CMatrix(MatrixSize, MatrixSize);
  HamiltonianMatrix->InitializeAll(0.);


  ///The main stuff: construct the Hamiltoninan matrix.

  myPrinter.Print(1,"Constructing Hamiltonian matrix.\n");
  

  ParametrizedCurve * myCurve = myConfiguration.GetKCurve();
  Potential * myPotential = myConfiguration.GetPotential();

  myPrinter.Print(3, "Potential: minX %13.10e maxX %13.10e\n", myPotential->GetMinX(), myPotential->GetMaxX());

  MultiTasker<WorkerData, void*> * myMultiTasker = NULL;
  ///Generate the Hamiltonian in parallell!
  if(numberOfParticles == 1)
	{
	  myMultiTasker = new MultiTasker<WorkerData, void*>(EvaluateSubMatrixOneParticle, myConfiguration.GetNumberOfThreads());
	}
  else if(numberOfParticles == 2)
	{
	  myMultiTasker = new MultiTasker<WorkerData, void*>(EvaluateSubMatrixTwoParticles, myConfiguration.GetNumberOfThreads());
	}
  else
	{
	  throw RLException("Improbable exception, this should never happen.");
	}
  if(myMultiTasker == NULL)
	{
	  throw RLException("Equally improbable.");
	}
	

  myMultiTasker->RegisterListener(&myPrinter);

  for(uint i = 0; i<MatrixSize; ++i)
	{
	  myMultiTasker->AddInput(WorkerData(HamiltonianMatrix,
										 myCurve,
										 myPotential,
										 myBasisFunctions,
										 numberOfGLPoints,
										 myConfiguration.GetSpecificUnits()->GetHbarTimesLambda(),
										 myConfiguration.GetSpecificUnits()->GetMassOverLambda2(),
										 0,
										 i,
										 i+1,
										 myConfiguration.GetCouplingCoefficient(),
										 MatrixSize)
							  );
	}
  myMultiTasker->LaunchThreads();
  myMultiTasker->PauseUntilOutputIsGenerated();
  myMultiTasker->DestroyThreads();
  delete myMultiTasker; myMultiTasker = NULL;

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

void * EvaluateSubMatrixOneParticle(WorkerData w)
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


void * EvaluateSubMatrixTwoParticles(WorkerData w)
{
  CMatrix * HamiltonianMatrix = w.HamiltonianMatrix;
  ParametrizedCurve * myCurve = w.myCurve;
  Potential * myPotential = w.myPotential;
  vector<BasisFunction> * myBasisFunctions = &w.myBasisFunctions;
  uint numberOfGLPoints = w.numberOfGLPoints;
  double hbarTimesLambda = w.hbarTimesLambda;
  double massOverLambda2 = w.massOverLambda2;
  uint m1 = w.m1, m2 = w.m2, n1 = w.n1, n2 = w.n2;
  double couplingCoefficient = w.couplingCoefficient;
  uint N = myBasisFunctions->size() * numberOfGLPoints; ///Number of points in each subspace.



  for(uint i = m1; i<m2; ++i)
	{
	  uint a = i/N;
	  uint b = i % N;

	  uint curvePointerA = a % numberOfGLPoints;
	  uint basisPointerA = a / numberOfGLPoints;
	  uint curveSegmentA = myCurve->SegmentIndexFromGLNumber(curvePointerA);
	  
	  ComplexDouble kA = myCurve->GetRuleValue(curveSegmentA, curvePointerA);
	  ComplexDouble wA = myCurve->GetRuleWeight(curveSegmentA, curvePointerA);


	  uint curvePointerB = b % numberOfGLPoints;
	  uint basisPointerB = b / numberOfGLPoints;
	  uint curveSegmentB = myCurve->SegmentIndexFromGLNumber(curvePointerB);
	  
	  ComplexDouble kB = myCurve->GetRuleValue(curveSegmentB, curvePointerB);
	  ComplexDouble wB = myCurve->GetRuleWeight(curveSegmentB, curvePointerB);

	  
	  for(uint j = n1; j<n2; ++j)
		{
		  uint c = j / N;
		  uint d = j % N;

		  if(c != d && a != b) //To speed things up.
			continue;

		  uint curvePointerC = c % numberOfGLPoints;
		  uint basisPointerC = c / numberOfGLPoints;
		  uint curveSegmentC = myCurve->SegmentIndexFromGLNumber(curvePointerC);
		  
		  ComplexDouble kC = myCurve->GetRuleValue(curveSegmentC, curvePointerC);
		  ComplexDouble wC = myCurve->GetRuleWeight(curveSegmentC, curvePointerC);

		  uint curvePointerD = d % numberOfGLPoints;
		  uint basisPointerD = d / numberOfGLPoints;
		  uint curveSegmentD = myCurve->SegmentIndexFromGLNumber(curvePointerD);
		  
		  ComplexDouble kD = myCurve->GetRuleValue(curveSegmentD, curvePointerD);
		  ComplexDouble wD = myCurve->GetRuleWeight(curveSegmentD, curvePointerD);

		  if(b==d) //2nd term
			{
			  HamiltonianMatrix->Element(i, j) += 
				sqrt(wA*wC) * myPotential->BasisIntegrate((*myBasisFunctions)[basisPointerA], 
														   (*myBasisFunctions)[basisPointerC], 
														   kA, kC);
			}

		  if(a==c) //3rd term
			{
			  HamiltonianMatrix->Element(i, j) += 
				sqrt(wB*wD) * myPotential->BasisIntegrate((*myBasisFunctions)[basisPointerB], 
														   (*myBasisFunctions)[basisPointerD], 
														   kB, kD);
			}

		  if(a==c && b==d) //1st term
			{
			  HamiltonianMatrix->Element(i, j) += 
				pow(hbarTimesLambda,2)/(2.*massOverLambda2) * (pow(kA,2) + pow(kB,2) + couplingCoefficient);
			}
		}
   }
   return NULL;
}

