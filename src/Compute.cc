#include "Compute.hh"



ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2)
{
  return myPotential->Evaluate(x) * exp(ComplexDouble(0,1)*(k1-k2));
}




int main(int argc, char *argv[])
{
  CommandLineInterpreter * myInterpreter = InitInterpreter(); 
  try
	{
	  myInterpreter->Initialize(argc, argv);
	  if(!myInterpreter->ReadFlaggedCommand("help").empty())
		{
		  myInterpreter->PrintHelp();
		  return 0;
		}
	}
  catch(CommandLineException e)
	{
	  myInterpreter->PrintHelp();
	  return 1;
	}

  int verbosityLevel = 0;
  if(myInterpreter->ReadFlaggedCommand("verbose").size()==1)
    verbosityLevel = atoi((myInterpreter->ReadFlaggedCommand("verbose").front()).c_str());

  unsigned int threads = 10;
  if(!myInterpreter->ReadFlaggedCommand("threads").empty())
      threads = atoi(myInterpreter->ReadFlaggedCommand("threads").front().c_str());



  int nKValues = 100;
  myPotential = new Potential();
  KPoints myKPoints(nKValues, 50, 2, 2);
  int xPoints = 2000;

  CMatrix HamiltonianMatrix(2*nKValues+1,2*nKValues+1);
  HamiltonianMatrix.InitializeAll(0.);
  
  foru(i, (int)myKPoints.GetPoints()->size())
	{
	  foru(j, (int)myKPoints.GetPoints()->size())
		{

		  HamiltonianMatrix.Element(i, j) += Integrator::Integrate( IntegrandValue, xPoints, myPotential->GetMinX(), myPotential->GetMaxX(), myKPoints.GetPoint(i), myKPoints.GetPoint(j));
		  ///Fix the discretization step stuff!!!
		  HamiltonianMatrix.Element(i, j)/=2*M_PI;
		  if ( i==j )
			{
			  HamiltonianMatrix.Element(i, j) += (pow(HBAR, 2)/(2.*MASS) * pow(myKPoints.GetPoint(i), 2))/(real(myKPoints.GetDeltaK(i)));
			}
		}
	}
  HamiltonianMatrix.MultiplyBy(myKPoints.GetStencilDeltaK());

  delete myPotential;
    if ( ! HamiltonianMatrix.IsHermitian(true) )
	  {
		throw RLException("The matrix was found to be non-hermitian.");
	  }
	
  EigenInformation myInfo = EigenvalueSolver::Solve(&HamiltonianMatrix);

  return 0;
}



CommandLineInterpreter * InitInterpreter()
{
  CommandLineInterpreter * myInterpreter = new CommandLineInterpreter();

  myInterpreter->AddCommandLineArgument(CommandLineArgument("help",0,false, "Displays a help message and quits."));

  myInterpreter->SetDescription("Compute eigenvalues for particle states in a potential.");
  return myInterpreter;
}

