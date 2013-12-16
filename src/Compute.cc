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


  OutputProcessor myProcessor(&myConfiguration);
  myProcessor.RegisterListener(&myPrinter);


  vector<BasisFunction> myBasisFunctions = myConfiguration.GetBasisFunctions();
  myPrinter.Print(5, "Using %d basis functions:", myBasisFunctions.size());
  for(vector<BasisFunction>::const_iterator it = myBasisFunctions.begin(); it!=myBasisFunctions.end(); ++it)
	{
	  myPrinter.Print(5, "%s %s", ((it==myBasisFunctions.begin())?"":","),it->GetName());
	}
  myPrinter.Print(5, ".\n");

  unsigned int numberOfGLPoints = myConfiguration.GetKCurve()->GetTotalNumberOfGLPoints();


  unsigned int numberOfBasisFunctions =  myBasisFunctions.size();
  unsigned int MatrixSize = numberOfGLPoints * numberOfBasisFunctions;
  myPrinter.Print(5, "Total number of GL points in k-space: %d, Matrix size: %d.\n", numberOfGLPoints, MatrixSize);


  CMatrix HamiltonianMatrix(MatrixSize, MatrixSize);
  HamiltonianMatrix.InitializeAll(0.);


  ///The main stuff: construct the Hamiltoninan matrix.

  myPrinter.Print(1,"Constructing Hamiltonian matrix.\n");
  

  ParametrizedCurve * myCurve = myConfiguration.GetKCurve();
  Potential * myPotential = myConfiguration.GetPotential();

  myPrinter.Print(3, "Potential: minX %13.10e maxX %13.10e\n", myPotential->GetMinX(), myPotential->GetMaxX());


  ///Generate the Hamiltonian in parallell!
  MultiTasker<WorkerData, void*> myMultiTasker(EvaluateSubMatrix, myConfiguration.GetNumberOfThreads());

  myMultiTasker.RegisterListener(&myPrinter);

  for(unsigned int i = 0; i<MatrixSize; ++i)
	{
	  myMultiTasker.AddInput(WorkerData(&HamiltonianMatrix,
										myCurve,
										myPotential,
										myBasisFunctions,
										numberOfGLPoints,
										myConfiguration.GetSpecificUnits()->GetHbarTimesLambda(),
										myConfiguration.GetSpecificUnits()->GetMassOverLambda2(),
										i,
										i+1,
										0,
										MatrixSize)
						);
	}
  myMultiTasker.LaunchThreads();
  myMultiTasker.PauseUntilOutputIsGenerated();
  myMultiTasker.DestroyThreads();
  
  ///Now we should have a Hamiltonian.


  ///Validate some basic properties of the Hamilton matrix, if they are expected.
  if(myConfiguration.GetExpectedMatrixType() == SymmetricMatrix)
	{
	  myPrinter.Print(1, "Validating symmetricity of matrix.\n");
	  if ( ! HamiltonianMatrix.IsSymmetric(true) )
		{
		  throw RLException("The matrix was found to be non-symmetric.");
		}
	}
  if(myConfiguration.GetExpectedMatrixType() == HermitianMatrix)
	{
	  myPrinter.Print(1, "Validating hermiticity of matrix.\n");
	  if( ! HamiltonianMatrix.IsHermitian(true))
		{
		  throw RLException("The matrix was found to be non-hermitian.");
		}
	}

  myPrinter.Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");
  EigenInformation myInfo = EigenvalueSolver::Solve(&HamiltonianMatrix);
  
  
  myPrinter.Print(2, "Saving results to files.\n");

  myProcessor.SetEigenInformation(&myInfo);
  myProcessor.WriteOutput();


  myPrinter.Print(2, "Launching plotters.\n");
  
  ExternalLauncher myLauncher(&myConfiguration);
  myLauncher.Launch();

  myPrinter.Print(2, "Done, exiting.\n");
  return 0;
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

void * EvaluateSubMatrix(WorkerData w)
{
  CMatrix * HamiltonianMatrix = w.HamiltonianMatrix;
  ParametrizedCurve * myCurve = w.myCurve;
  Potential * myPotential = w.myPotential;
  vector<BasisFunction> * myBasisFunctions = &w.myBasisFunctions;
  unsigned int numberOfGLPoints = w.numberOfGLPoints;
  double hbarTimesLambda = w.hbarTimesLambda;
  double massOverLambda2 = w.massOverLambda2;
  unsigned int m1 = w.m1, m2 = w.m2, n1 = w.n1, n2 = w.n2;

  for(unsigned int i = m1; i<m2; ++i)
	{
	  unsigned int curvePointerA = i % numberOfGLPoints;
	  unsigned int basisPointerA = i / numberOfGLPoints;
	  unsigned int curveSegmentA = myCurve->SegmentIndexFromGLNumber(curvePointerA);
	  
	  ComplexDouble kA = myCurve->GetRuleValue(curveSegmentA, curvePointerA);
	  ComplexDouble wA = myCurve->GetRuleWeight(curveSegmentA, curvePointerA);
	  
	  //ComplexDouble preFactorA = myBasisFunctions[basisPointerA].GetPreFactor();
	  
	  for(unsigned int j = n1; j<n2; ++j)
		{
		  unsigned int curvePointerB = j % numberOfGLPoints;
		  unsigned int basisPointerB = j / numberOfGLPoints;
		  unsigned int curveSegmentB = myCurve->SegmentIndexFromGLNumber(curvePointerB);
		  
		  ComplexDouble kB = myCurve->GetRuleValue(curveSegmentB, curvePointerB);
		  ComplexDouble wB = myCurve->GetRuleWeight(curveSegmentB, curvePointerB);
		  
		  //ComplexDouble preFactorB = myBasisFunctions[basisPointerB].GetPreFactor();
		  
		  HamiltonianMatrix->Element(i, j) += sqrt(wA*wB)*
			myPotential->BasisIntegrate((*myBasisFunctions)[basisPointerA], 
										(*myBasisFunctions)[basisPointerB], 
										kA, kB) ;
		}
	HamiltonianMatrix->Element(i,i) += pow(hbarTimesLambda,2)/(2.*massOverLambda2) * pow(kA, 2);
   }
   return NULL;
}
