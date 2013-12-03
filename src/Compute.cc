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

  PrintPotentialToFile(myConfiguration.GetPotentialFile(), myConfiguration.GetPotential());
  PrintPotentialPrecisionToFile(myConfiguration.GetPotentialPrecisionFile(), myConfiguration.GetPotential());
  PrintParametrizedCurveToFile(myConfiguration.GetKCurveFile(), myConfiguration.GetKCurve());
  PrintKFoundToFile(myConfiguration.GetKFoundFile(), &myInfo);



  myPrinter.Print(2, "Launching plotters.\n");
  if(myConfiguration.GetAutoPlotPotential())
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./PotentialPlot.sh \"%s\" \"%s\"", myConfiguration.GetPotentialFile(), myConfiguration.GetPotentialPrecisionFile());
	  system(buffer);
	}
  if(myConfiguration.GetAutoPlotKCurve())
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./KPlot.sh \"%s\" \"%s\"", myConfiguration.GetKCurveFile(), myConfiguration.GetKFoundFile());
	  system(buffer);
	}

  PrintInterestingKPoints(&myInfo,myConfiguration.GetKCurve());


  myPrinter.Print(2, "Done, exiting.\n");
  return 0;
}


void PrintPotentialToFile(const char * fileName, const Potential * potential)
{
  FILE * fout = fopen(fileName, "w");

  double potentialLength = potential->GetMaxX() - potential->GetMinX();
  if(potentialLength < 0)
	throw RLException("The potential length was calculated to be < 0. Something is wrong.");
  
  double scaleFactor = 0.2; ///How much zero spacing to add at the end of the potential.

  vector<pair<double, double> >  plottingPoints = potential->GetPlottingPoints();

  fprintf(fout, "#Potential function\n");
  fprintf(fout, "%+13.10e %+13.10e\n", potential->GetMinX() - scaleFactor*potentialLength, 0.);
  for(vector<pair<double, double> >::const_iterator it = plottingPoints.begin(); it!=plottingPoints.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", it->first, it->second);
	}
  fprintf(fout, "%+13.10e %+13.10e\n", potential->GetMaxX() + scaleFactor*potentialLength,0.0);
 
  fclose(fout);
  fout = NULL; 
}

void PrintPotentialPrecisionToFile(const char * fileName, const Potential * potential)
{
  FILE * fout = fopen(fileName, "w");

  vector<pair<double, double> >  plottingPoints = potential->GetPrecisionPoints();
  fprintf(fout, "#Potential precision\n");
  for(vector<pair<double, double> >::const_iterator it = plottingPoints.begin(); it!=plottingPoints.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", it->first, it->second);
	}
 
  fclose(fout);
  fout = NULL; 
}

void PrintParametrizedCurveToFile(const char * fileName, const ParametrizedCurve * toPrint)
{
  vector<ComplexDouble> printVector;
  for(unsigned int i = 0; i<toPrint->GetNumberOfSegments(); ++i)
	{
	  const vector<pair<ComplexDouble, ComplexDouble> > * rule = toPrint->GetSegmentRule(i);
	  for(vector<pair<ComplexDouble, ComplexDouble> >::const_iterator it = rule->begin(); it != rule->end(); ++it)
		{
		  printVector.push_back(it->first);
		}
	}
  PrintKCurveToFile(fileName, printVector);
}

void PrintInterestingKPoints(const EigenInformation * toPrint, const ParametrizedCurve * filter)
{
  vector<ComplexDouble> filterVector;
  for(unsigned int i = 0; i<filter->GetNumberOfSegments(); ++i)
	{
	  const vector<pair<ComplexDouble, ComplexDouble> > * rule = filter->GetSegmentRule(i);
	  for(vector<pair<ComplexDouble, ComplexDouble> >::const_iterator it = rule->begin(); it != rule->end(); ++it)
		{
		  filterVector.push_back(it->first);
		}
	}

  double minDistance = 1E99;
  double maxDistance = -1;
  for(unsigned int i = 1; i<filterVector.size(); ++i)
	{
	  double localDist = abs(filterVector[i-1] - filterVector[i]);
	  minDistance = MIN(minDistance, localDist);
	  maxDistance = MAX(maxDistance, localDist);
	}

  vector<ComplexDouble> printVector;
  for(vector<ComplexDouble>::const_iterator it = toPrint->Eigenvalues.begin(); it!=toPrint->Eigenvalues.end(); ++it)
	{
	  ComplexDouble kToPrint = sqrt((*it)*(double)2.*(double)MASSOVERC2)/HBARC;

	  ///If numerical stability is mean to us, then rotate.
	  if( (abs(imag(kToPrint)) > 1E2*abs(real(kToPrint))  && imag(kToPrint) < 0))
		{
		  kToPrint *= -1.0;
		}

	  ///Now apply filter rule.
	  if(imag(kToPrint) > minDistance )
		printVector.push_back(kToPrint);
	  else if(real(kToPrint) > minDistance)
		{
		  bool tooClose = false;
		  for(unsigned int i = 0; i<filterVector.size(); ++i)
			{
			  double d1 = abs(filterVector[i] - kToPrint);
			  double d2 = 0;
			  if(i>0)
				d2 += abs(filterVector[i-1]-filterVector[i]);
			  if(i+1<filterVector.size())
				d2 += abs(filterVector[i+1]-filterVector[i]);
			  if(i > 0 && i+1 < filterVector.size())
				d2 /= 2;

			  d2 *= 0.6;

			  if(d1 < d2)
				{
				  tooClose = true;
				  break;
				}
			}
		  if(!tooClose)
			printVector.push_back(kToPrint);
		}
	}

  for(vector<ComplexDouble>::const_iterator it = printVector.begin(); it!=printVector.end(); ++it)
	{
	  if(imag(*it) > minDistance && abs(arg(*it)-PI/2) < 1E-2 )
		{
		  printf("Bound state: %+6.10lfi\n", imag(*it));
		}
	}
  for(vector<ComplexDouble>::const_iterator it = printVector.begin(); it!=printVector.end(); ++it)
	{
	  if(!(imag(*it) > minDistance && abs(arg(*it)-PI/2) < 1E-2 ))
		{
		  printf("Resonant state: %+6.10lf %+6.10lfi\n", real(*it), imag(*it));
		}
	}
 
}

void PrintKFoundToFile(const char * fileName, const EigenInformation * toPrint)
{
  vector<ComplexDouble> printVector;
  for(vector<ComplexDouble>::const_iterator it = toPrint->Eigenvalues.begin(); it!=toPrint->Eigenvalues.end(); ++it)
	{
	  ComplexDouble kToPrint = sqrt((*it)*(double)2.*(double)MASSOVERC2)/HBARC;

	  ///If numerical stability is mean to us, then rotate.
	  if( (abs(imag(kToPrint)) > 1E2*abs(real(kToPrint))  && imag(kToPrint) < 0))
		{
		  kToPrint *= -1.0;
		}

	  printVector.push_back(kToPrint);
	}
  
  PrintKCurveToFile(fileName, printVector);
}

void PrintKCurveToFile(const char * fileName, const vector<ComplexDouble> & toPrint)
{
  FILE * fout = fopen(fileName, "w");
  fprintf(fout, "#KCurve\n");
  for(vector<ComplexDouble>::const_iterator it = toPrint.begin(); it!=toPrint.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", real(*it), imag(*it));
	}
  fclose(fout);
  fout = NULL;
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
		  
		  HamiltonianMatrix->Element(i, j) += ComplexDouble(1./(2.*PI),0)*
			sqrt(wA*wB)*
			myPotential->BasisIntegrate((*myBasisFunctions)[basisPointerA], 
										(*myBasisFunctions)[basisPointerB], 
										kA, kB) ;
		}
	HamiltonianMatrix->Element(i,i) += pow(HBARC,2)/(2.*MASSOVERC2) * pow(kA, 2);
   }
   return NULL;
}
