#include "Compute.hh"


/*
ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2)
{
  return myPotential->Evaluate(x) * exp(ComplexDouble(0,1)*(k1-k2));
}
*/


ComplexDouble ExpOfInnerProduct(ComplexDouble k1, ComplexDouble k2)
{
  return k1 + k2;
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

  Potential myPotential;

  double kPMax = 2*kValues + 1;

  myPrinter->Print(2, "Computing Legendre rule.\n");
  vector<pair<double, double> > myLegendreRule = LegendreRule::GetRule(kPMax);



  CMatrix HamiltonianMatrix(kPMax,kPMax);
  HamiltonianMatrix.InitializeAll(0.);


  myPrinter->Print(1,"Constructing Hamiltonian matrix.\n");
  for(unsigned int i = 0; i<kPMax; ++i)
	{
	  ComplexDouble ki = kCurve.Evaluate(myLegendreRule[i].first);
	  double wi = myLegendreRule[i].second;
	  for(unsigned int j = 0; j<kPMax; ++j)
		{
		  ComplexDouble kj = kCurve.Evaluate(myLegendreRule[j].first);
		  double wj = myLegendreRule[j].second;

		  HamiltonianMatrix.Element(i, j) += sqrt(wi*wj)*
			myPotential.FastExpIntegrate(ExpOfInnerProduct(kj, ki));
		}
	  HamiltonianMatrix.Element(i,i) -= pow(HBAR,2.)/(2.*MASS) * ki*ki;
	}

  myPrinter->Print(1, "Validating symmetricity of matrix.\n");
  
  if ( ! HamiltonianMatrix.IsSymmetric(true) )
	{
	  throw RLException("The matrix was found to be non-symmetric.");
	}
  
  myPrinter->Print(1, "Solving for eigenvalues and eigenvectors of the Hamiltonian.\n");

  EigenInformation myInfo = EigenvalueSolver::Solve(&HamiltonianMatrix);


  myPrinter->Print(2, "Cleaning up.\n");
  delete myPrinter;

  return 0;
}



CommandLineInterpreter * InitInterpreter()
{
  CommandLineInterpreter * myInterpreter = new CommandLineInterpreter();
  TMP_LIST_STR(kCutoffDescription,"double");
  TMP_LIST_STR(defaultCutoff, "50");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kCutoff", 1, false, "Cutoff for the k-values.", kCutoffDescription, defaultCutoff));

  TMP_LIST_STR(kValuesDescription, "double");
  TMP_LIST_STR(kValuesDefault, "100");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kValues", 1, false, "Number of k-values to use.", kValuesDescription, kValuesDefault));

  TMP_LIST_STR(kMidDescription, "double");
  TMP_LIST_STR(kMidDefault, "5");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kMid", 1, false, "Midpoint of dent in k-plane curve.", kMidDescription, kMidDefault));

  TMP_LIST_STR(kDepthDescription, "double");
  TMP_LIST_STR(kDepthDefault, "3");
  myInterpreter->AddCommandLineArgument(CommandLineArgument("kDepth", 1, false, "Depth of dent in k-plane curve.", kDepthDescription, kDepthDefault));

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

