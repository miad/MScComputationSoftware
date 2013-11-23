#include "Compute.hh"

/*
ComplexDouble ExpOfInnerProduct(ComplexDouble k1, ComplexDouble k2)
{
  if(invertInnerProduct)
	return (k1 - k2);

  return (k1 + k2);
}
*/

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

  string configFile = myInterpreter->ReadFlaggedCommandStrict("configFile").front().c_str();
  ComputeConfig myConfiguration;
  if( access(configFile.c_str(), F_OK) == -1 )
	{
	  myConfiguration.WriteFile(configFile.c_str());
	}
  else
	{
	  myConfiguration.ReadFile(configFile.c_str());
	}
  
  VerbosePrinter * myPrinter = new VerbosePrinter(myConfiguration.GetVerbosityLevel());
  
  myPrinter->Print(3, "Initializing\n");

  vector<BasisFunction> myBasisFunctions = myConfiguration.GetBasisFunctions();

  unsigned int numberOfGLPoints = myConfiguration.GetKCurve()->GetTotalNumberOfGLPoints();
  unsigned int numberOfBasisFunctions =  myBasisFunctions.size();
  unsigned int MatrixSize = numberOfGLPoints * numberOfBasisFunctions;



  CMatrix HamiltonianMatrix(MatrixSize, MatrixSize);
  HamiltonianMatrix.InitializeAll(0.);

  myPrinter->Print(1,"Constructing Hamiltonian matrix.\n");

  unsigned int hFactor = MAX(numberOfBasisFunctions, numberOfGLPoints);

  ParametrizedCurve * myCurve = myConfiguration.GetKCurve();
  Potential * myPotential = myConfiguration.GetPotential();

  for(unsigned int i = 0; i<MatrixSize; ++i)
	{
	  unsigned int curvePointerA = i % hFactor;
	  unsigned int basisPointerA = i / hFactor;
	  unsigned int curveSegmentA = myCurve->SegmentIndexFromGLNumber(curvePointerA);

	  ComplexDouble kA = myCurve->GetRuleValue(curveSegmentA, curvePointerA);
	  ComplexDouble wA = myCurve->GetRuleWeight(curveSegmentA, curvePointerA);
		
	  for(unsigned int j = 0; j<MatrixSize; ++j)
		{
		  unsigned int curvePointerB = j % hFactor;
		  unsigned int basisPointerB = j / hFactor;
		  unsigned int curveSegmentB = myCurve->SegmentIndexFromGLNumber(curvePointerB);
		  
		  ComplexDouble kB = myCurve->GetRuleValue(curveSegmentB, curvePointerB);
		  ComplexDouble wB = myCurve->GetRuleWeight(curveSegmentB, curvePointerB);
		  
		  HamiltonianMatrix.Element(i, j) += ComplexDouble(1./(2.*PI),0)*
			sqrt(wA*wB)
			*myPotential->BasisIntegrate(myBasisFunctions
		}
	  HamiltonianMatrix.Element(i,i) += pow(HBARC,2)/(2.*MASSOVERC2) * pow(kA, 2);
	}


  //myPrinter->Print(12,HamiltonianMatrix.ToString().c_str());


  if(myConfiguration.GetExpectedMatrixType() == SymmetricMatrix)
	{
	  myPrinter->Print(1, "Validating symmetricity of matrix.\n");
	  if ( ! HamiltonianMatrix.IsSymmetric(true) )
		{
		  throw RLException("The matrix was found to be non-symmetric.");
		}
	}
  if(myConfiguration.GetExpectedMatrixType() == HermitianMatrix)
	{
	  myPrinter->Print(1, "Validating hermiticity of matrix.\n");
	  if( ! HamiltonianMatrix.IsHermitian(true))
		{
		  throw RLException("The matrix was found to be non-hermitian.");
		}
	}

  
  myPrinter->Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");

  EigenInformation myInfo = EigenvalueSolver::Solve(&HamiltonianMatrix);

  /*
  for(vector<ComplexDouble>::const_iterator it = myInfo.Eigenvalues.begin(); it!=myInfo.Eigenvalues.end(); ++it)
	{
	  cout << *it << endl;
	}
  */

  PrintDataToFile(myPrinter, dataFile, myInfo, kValuesOnCurve, myConfiguration.GetPotential().GetPotentialPoints());


  myPrinter->Print(2, "Cleaning up.\n");
  delete myPrinter;

  return 0;
}

void PrintDataToFile(VerbosePrinter * myPrinter, const string fileName, const EigenInformation & data, const vector<ComplexDouble> & kValuesOnCurve, const list<Interval> & potentialIntervals)
{
  FILE * fout = fopen(fileName.c_str(), "w");
  fprintf(fout, "#This file was created by 'Compute', a computation program, for plotting stuff related to Rikard Lundmark's Master thesis.\n");

  ///Print a couple of datasets, to be able to visualize with GNUPLOT


  fprintf(fout, "\"k-values in basis\"\n");
  for(vector<ComplexDouble>::const_iterator it = kValuesOnCurve.begin(); it!=kValuesOnCurve.end(); ++it)
	{
	  //if(real(*it)>= 0)
		fprintf(fout, "%13.6e %13.6e\n",real(*it),imag(*it));
	}
  fprintf(fout, "\n\n\n");



  fprintf(fout, "\"Potential\"\n");

  double lastPoint = -10;
  for(list<Interval>::const_iterator it = potentialIntervals.begin(); it!=potentialIntervals.end(); ++it)
	{
	  if(!DBL_EQUAL(it->x1,lastPoint))
		{
		  fprintf(fout, "%13.6e %13.6e\n", lastPoint+EPS, 0.);
		  fprintf(fout, "%13.6e %13.6e\n", it->x1-EPS,0.);
		}
	  fprintf(fout, "%13.6e %13.6e\n", it->x1, it->y);
	  fprintf(fout, "%13.6e %13.6e\n", it->x2, it->y);
	  lastPoint = it->x2;
	}
  fprintf(fout, "%13.6e %13.6e\n", potentialIntervals.back().x2, 0.);
  fprintf(fout, "%13.6e %13.6e\n", 10., 0.);


  fprintf(fout, "\n\n\n");
  
  fprintf(fout, "\"k-values\"\n");
  for(vector<ComplexDouble>::const_iterator it = data.Eigenvalues.begin(); it!=data.Eigenvalues.end(); ++it)
	{
	  ComplexDouble kToPrint = sqrt((*it)*(double)2.*(double)MASSOVERC2)/HBARC;
	  ///Transform to uhp if numerical stability screwed us over...
	  if( (abs(imag(kToPrint)) > 1E1*abs(real(kToPrint))  && imag(kToPrint) < 0))
		{
		  if(myPrinter != NULL)
			{
			  myPrinter->Print(2,"Artificial rotation of a point to uhp.\n");
			}
		  kToPrint *= ComplexDouble(-1,0);
		}


	  fprintf(fout, "%13.6e, %13.6e\n", real(kToPrint), imag(kToPrint));
	}
  //  printf("%13.6e + %13.6ei    %13.6e + %13.6ei    %13.6e + %13.6ei\n",real(sqrt(ComplexDouble(-1,0))),imag(sqrt(ComplexDouble(-1,0))), real(sqrt(ComplexDouble(0,1))), imag(sqrt(ComplexDouble(0,1))), real(sqrt(ComplexDouble(0,-1))), imag(sqrt(ComplexDouble(0,-1))));


  fclose(fout);
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

