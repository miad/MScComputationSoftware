#include "Compute.hh"


ComplexDouble ExpOfInnerProduct(ComplexDouble k1, ComplexDouble k2)
{
  if(invertInnerProduct)
	return (k1 - k2);

  return (k1 + k2);
}

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

  ///Read out parameters from the command line interpreter (or default ones, if they were not specified.
  int verbosityLevel = atoi((myInterpreter->ReadFlaggedCommandStrict("verbose").front()).c_str());
  unsigned int threads = atoi(myInterpreter->ReadFlaggedCommandStrict("threads").front().c_str());
  double kCutoff = atof(myInterpreter->ReadFlaggedCommandStrict("kCutoff").front().c_str());
  double kMid = atof(myInterpreter->ReadFlaggedCommandStrict("kMid").front().c_str());
  double kDepth = atof(myInterpreter->ReadFlaggedCommandStrict("kDepth").front().c_str());
  unsigned int kValues = atof(myInterpreter->ReadFlaggedCommandStrict("kValues").front().c_str());
  string dataFile = myInterpreter->ReadFlaggedCommandStrict("dataFile").front().c_str();
  string potentialFile = myInterpreter->ReadFlaggedCommandStrict("potentialFile").front().c_str();

  invertInnerProduct = !myInterpreter->ReadFlaggedCommand("flipInner").empty();

  VerbosePrinter * myPrinter = new VerbosePrinter(verbosityLevel);

  myPrinter->Print(3, "Initializing k-curve.\n");
  ///Create the standard Berggren curve.
  ParametrizedCurve kCurve;
  kCurve.AddValue(-kCutoff);
  kCurve.AddValue(-2*kMid);
  kCurve.AddValue(-kMid + (kDepth*ComplexDouble(0, 1)));
  kCurve.AddValue(kMid  - (kDepth*ComplexDouble(0,1)));
  kCurve.AddValue(2*kMid);
  kCurve.AddValue(kCutoff);

  Potential myPotential(potentialFile);

  unsigned int kPMax = kValues*kCurve.GetNumberOfSegments();

  myPrinter->Print(2, "Computing Legendre rule.\n");
  vector<pair<double, double> > myLegendreRule = LegendreRule::GetRule(kValues);



  CMatrix HamiltonianMatrix(kPMax,kPMax);
  HamiltonianMatrix.InitializeAll(0.);

  vector<ComplexDouble> kValuesOnCurve;
  for(unsigned int j = 0; j<kCurve.GetNumberOfSegments(); ++j)
	{
	  for(unsigned int i = 0; i<kValues; ++i)
		{
		  kValuesOnCurve.push_back(kCurve.SegmentEvaluate(j, myLegendreRule[i].first));
		}
	}


  myPrinter->Print(1,"Constructing Hamiltonian matrix.\n");
  for(unsigned int i = 0; i<kPMax; ++i)
	{
	  ComplexDouble ki = kValuesOnCurve[i];
	  double wi = myLegendreRule[i%kValues].second;
	  for(unsigned int j = 0; j<kPMax; ++j)
		{
		  ComplexDouble kj = kValuesOnCurve[j];
		  double wj = myLegendreRule[j%kValues].second;

		  HamiltonianMatrix.Element(i, j) += ComplexDouble(1./(2.*PI*HBAR),0)*
			sqrt(wi*wj)*
			sqrt(kCurve.GetSegmentDerivative(i/kValues)*kCurve.GetSegmentDerivative(j/kValues))*
			myPotential.FastExpIntegrate(ExpOfInnerProduct(kj, ki)/ComplexDouble(HBAR,0));
		}
	  HamiltonianMatrix.Element(i,i) += 1./(2.*MASS) * ki*ki;
	}

  if(verbosityLevel > 10)
	{
	  myPrinter->Print(12, "Hamiltonian matrix: \n");
	  myPrinter->Print(12,HamiltonianMatrix.ToString().c_str());
	}


  if(!invertInnerProduct)
	{
	  myPrinter->Print(1, "Validating symmetricity of matrix.\n");
	  if ( ! HamiltonianMatrix.IsSymmetric(true) )
		{
		  throw RLException("The matrix was found to be non-symmetric.");
		}
	}
  else
	{
	  myPrinter->Print(1, "Inner product inverted: not performing symmetry check.\n");
	}
  
  myPrinter->Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");

  EigenInformation myInfo = EigenvalueSolver::Solve(&HamiltonianMatrix);

  /*
  for(vector<ComplexDouble>::const_iterator it = myInfo.Eigenvalues.begin(); it!=myInfo.Eigenvalues.end(); ++it)
	{
	  cout << *it << endl;
	}
  */

  PrintDataToFile(dataFile, myInfo, kValuesOnCurve, myPotential.GetPotentialPoints());


  myPrinter->Print(2, "Cleaning up.\n");
  delete myPrinter;

  return 0;
}

void PrintDataToFile(const string fileName, const EigenInformation & data, const vector<ComplexDouble> & kValuesOnCurve, const list<Interval> & potentialIntervals)
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
	  ComplexDouble kToPrint = sqrt((*it)*(double)2.*(double)MASS);
	  ///Transform to uhp if numerical stability screwed us over...
	  if( (abs(imag(kToPrint)) > 1E1*abs(real(kToPrint))  && imag(kToPrint) < 0))
		kToPrint *= ComplexDouble(-1,0);


	  fprintf(fout, "%13.6e, %13.6e\n", real(kToPrint), imag(kToPrint));
	}
  //  printf("%13.6e + %13.6ei    %13.6e + %13.6ei    %13.6e + %13.6ei\n",real(sqrt(ComplexDouble(-1,0))),imag(sqrt(ComplexDouble(-1,0))), real(sqrt(ComplexDouble(0,1))), imag(sqrt(ComplexDouble(0,1))), real(sqrt(ComplexDouble(0,-1))), imag(sqrt(ComplexDouble(0,-1))));


  fclose(fout);
}



CommandLineInterpreter * InitInterpreter()
{
  CommandLineInterpreter * myInterpreter = new CommandLineInterpreter();
  TMP_LIST_STR(kCutoffDescription,"double");
  TMP_LIST_STR(defaultCutoff, "50");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kCutoff", 1, false, "Cutoff for the k-values.", kCutoffDescription, defaultCutoff));

  TMP_LIST_STR(kValuesDescription, "double");
  TMP_LIST_STR(kValuesDefault, "50");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kValues", 1, false, "Number of k-values per segment to use.", kValuesDescription, kValuesDefault));

  TMP_LIST_STR(kMidDescription, "double");
  TMP_LIST_STR(kMidDefault, "5");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kMid", 1, false, "Midpoint of dent in k-plane curve.", kMidDescription, kMidDefault));

  TMP_LIST_STR(kDepthDescription, "double");
  TMP_LIST_STR(kDepthDefault, "0.2");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kDepth", 1, false, "Depth of dent in k-plane curve.", kDepthDescription, kDepthDefault));

  TMP_LIST_STR(dataFileDescription, "file name");
  TMP_LIST_STR(dataFileDefault,"compute_output.dat");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("dataFile", 1, false, "Data file containing output data.", dataFileDescription, dataFileDefault));

  TMP_LIST_STR(potentialFileDescription, "file name");
  TMP_LIST_STR(potentialFileDefault,"potential.dat");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("potentialFile", 1, false, "Input file specifying the (piecewise) potential.", potentialFileDescription, potentialFileDefault));


  myInterpreter->AddCommandLineArgument(CommandLineArgument("flipInner", 0, false, "Use to flip sign of k2 in inner product (will cause non-symmetric matrix)."));


  TMP_LIST_STR(verboseDefault,"0");
  TMP_LIST_STR(verboseDescription,"positive integer");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("verbose", 1, false, "How verbose the program should be during the execution process.", verboseDescription, verboseDefault));



  TMP_LIST_STR(threadsDefault,"5");
  TMP_LIST_STR(threadsDescription,"positive integer");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("threads", 1, false, "How many threads to use in multi-threaded sections of the program.", threadsDescription, threadsDefault));

  myInterpreter->AddCommandLineArgument(CommandLineArgument("help",0,false, "Displays a help message and quits."));

  myInterpreter->SetDescription("Compute eigenvalues for particle states in a potential.");
  return myInterpreter;
}

