#include "ComputeConfig.hh"

ComputeConfig::ComputeConfig()
  :numberOfThreads(1), verbosityLevel(0), autoPlotPotential(false), autoPlotKCurve(false),
   minWavefunctionX(-10), maxWavefunctionX(10), wavefunctionStepsizeX(0.01)
{
  strcpy(kCurveFile, "KCurve.dat");
  strcpy(kFoundFile, "KFound.dat");
  strcpy(potentialFile,"Potential.dat");
  strcpy(potentialPrecisionFile,"PotentialPrecision.dat");
  strcpy(interestingPointsFile, "InterestingPoints.dat");
  strcpy(wavefunctionFile, "Wavefunctions.dat");
  PiecewiseConstantPotential * stdPotential = new PiecewiseConstantPotential();
  stdPotential->AddValue(-3, -1, 10);
  stdPotential->AddValue(-1, 1, -225);
  stdPotential->AddValue(1, 3, 10);
  stdPotential->RecomputeLegendreRules();

  potential = stdPotential;


  kCurve = new ParametrizedCurve(-1,1);
  kCurve->AddValue(0.0);
  kCurve->AddValue(ComplexDouble(3., -0.2));
  kCurve->AddValue(ComplexDouble(8.,0.));
  kCurve->AddValue(ComplexDouble(30.,0.));
  kCurve->AddGLPoints(25);
  kCurve->AddGLPoints(20);
  kCurve->AddGLPoints(30);
  kCurve->ComputeGaussLegendre();


  basisFunctions.push_back(BasisFunction("sqrt(2)*sin(k*x)"));
  basisFunctions.push_back(BasisFunction("sqrt(2)*cos(k*x)"));

}

ComputeConfig::ComputeConfig(const char * fileName)
{
  ReadFile(fileName);
}

ComputeConfig::~ComputeConfig()
{
  if(potential != NULL)
	{
	  delete potential;
	  potential = NULL;
	}
  if(kCurve != NULL)
	{
	  delete kCurve;
	  kCurve = NULL;
	}
  if(potential != NULL)
	{
	  delete potential;
	  potential = NULL;
	}
}

void ComputeConfig::ReadFile(const char * fileName)
{
  Config cfg;
  cfg.readFile(fileName);
  Setting & root = cfg.getRoot();

  delete potential;
  potential = NULL;

  kCurve->Clear();
  basisFunctions.clear();


  ///Validate version.
  double version;
  if( root.lookupValue("Version", version) )
	{
	  if( ! DBL_EQUAL(version, CONFIG_FILE_VERSION) )
		throw RLException("Invalid configuration file version: was %f, required %f", version, CONFIG_FILE_VERSION);
	}
  else
	{
	  throw RLException("Could not find property 'Version' in configuration file.");
	}

  if( ! root.exists("Program") || ! root["Program"].isGroup())
	{
	  throw RLException("The group 'Program' was not properly defined in the config file.");
	}

  Setting & program = root["Program"];
  if(!program.exists("AutoLaunch") || ! program["AutoLaunch"].isGroup())
	{
	  throw RLException("The group 'AutoLaunch' was not properly defined in the config file.");
	}

  if(! program.lookupValue("VerbosityLevel", verbosityLevel))
	{
	  throw RLException("The value 'VerbosityLevel' was not defined.");
	}

  if(! program.lookupValue("NumberOfThreads", numberOfThreads))
	{
	  throw RLException("The value 'NumberOfThreads' was not defined.");
	}

  Setting & autolaunch = program["AutoLaunch"];
  if(!autolaunch.exists("Gnuplot") || ! autolaunch["Gnuplot"].isGroup())
	{
	  throw RLException("The group 'Gnuplot' was not properly defined in the config file.");
	}

  Setting & gnuplot = autolaunch["Gnuplot"];

  if( !gnuplot.lookupValue("Potential", autoPlotPotential))
	{
	  throw RLException("The value 'Potential' in 'AutoLaunch' was not set appropriately.");
	}

  if( !gnuplot.lookupValue("KCurve", autoPlotKCurve))
	{
	  throw RLException("The value 'KCurve' in 'AutoLaunch' was not set appropriately.");
	}

  if( !root.exists("Output") || ! root["Output"].isGroup())
	{
	  throw RLException("The value 'Output' was not set appropriately.");
	}

  Setting & output = root["Output"];
  string temp;
  if( output.lookupValue("KCurve", temp) )
	{
	  strcpy(kCurveFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'KCurve' was not set appropriately.");
	}

  if( output.lookupValue("KFound", temp) )
	{
	  strcpy(kFoundFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'KFound' was not set appropriately.");
	}

  if( output.lookupValue("Potential", temp) )
	{
	  strcpy(potentialFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'Potential' was not set appropriately.");
	}

  if( output.lookupValue("PotentialPrecision", temp) )
	{
	  strcpy(potentialPrecisionFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'PotentialPrecision' was not set appropriately.");
	}

  if( output.lookupValue("InterestingPoints", temp) )
	{
	  strcpy(interestingPointsFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'InterestingPoints' was not set appropriately.");
	}

  if( output.lookupValue("Wavefunctions", temp) )
	{
	  strcpy(wavefunctionFile, temp.c_str());
	}
  else
	{
	  throw RLException("The value 'Wavefunctions' was not set appropriately.");
	}

  if( !output.exists("WavefunctionProperties") || !output["WavefunctionProperties"].isGroup() )
	{
	  throw RLException("WavefunctionProperties not properly defined.");
	}

  Setting & wfProp = output["WavefunctionProperties"];
  if( ! wfProp.lookupValue("MinX", minWavefunctionX))
	{
	  throw RLException("Could not find MinX in WavefunctionProperties.");
	}
  if( ! wfProp.lookupValue("MaxX", maxWavefunctionX))
	{
	  throw RLException("Could not find MaxX in WavefunctionProperties.");
	}
  if( ! wfProp.lookupValue("DeltaX", wavefunctionStepsizeX))
	{
	  throw RLException("Could not find DeltaX in WavefunctionProperties.");
	}



  if( !root.exists("Computation") || !root["Computation"].isGroup())
	{
	  throw RLException("'Computation' was not appropriately defined as a group.");
	}

  Setting & computation = root["Computation"];


  
  if( ! computation.lookupValue("ExpectedMatrixType", temp) )
	{
	  throw RLException("Expected matrix type not set properly.");
	}
  bool found = false;
  if(strcmp(temp.c_str(), "General") == 0) 
	{
	  matrixType = GeneralMatrix;
	  found = true;
	}
  if(strcmp(temp.c_str(), "Symmetric") == 0)
	{
	  matrixType = SymmetricMatrix;
	  found = true;
	}
  if(strcmp(temp.c_str(), "Hermitian") == 0 )
	{
	  matrixType = HermitianMatrix;
	  found = true;
	}
  if(!found)
	throw RLException("Invalid option for ExpectedMatrixType.");



  if( !computation.exists("BasisFunctions") || !computation["BasisFunctions"].isArray())
	{
	  throw RLException("Basis functions not properly defined.");
	}

  Setting & bf = computation["BasisFunctions"];
  for(int i = 0; i<bf.getLength(); ++i)
	{
	  basisFunctions.push_back(BasisFunction(bf[i]));
	}

  if( !computation.exists("Potential") || !computation["Potential"].isGroup())
	{
	  throw RLException("Potential not properly defined.");
	}

  Setting & poten = computation["Potential"];
  
  if( ! poten.lookupValue("Type", temp))
	{
	  throw RLException("Potential type not properly defined.");
	}



  if(strcmp(temp.c_str(), "PiecewiseConstant") == 0)
	{
	  PiecewiseConstantPotential * locPot = new PiecewiseConstantPotential();
	  int potPrec;
	  if( ! poten.lookupValue("Precision", potPrec))
		{
		  throw RLException("Could not find precision setting in settings file.");
		}
	  locPot->SetPrecision(potPrec);

	  if( !poten.exists("Values") || ! poten["Values"].isList())
		{
		  throw RLException("Could not find 'Values'.");
		}
	  
	  Setting & vval = poten["Values"];
	  for(int i = 0; i<vval.getLength(); ++i)
		{
		  double x1, x2, y;
		  if( ! vval[i].exists("Interval") || ! vval[i]["Interval"].isArray() || vval[i]["Interval"].getLength() != 2)
			{
			  throw RLException("Potential interval #%d was not set correctly.", i);
			}
		  x1 = vval[i]["Interval"][0];
		  x2 = vval[i]["Interval"][1];
		  
		  if( !vval[i].lookupValue("Value", y))
			{
			  throw RLException("Potential value in interval #%d was not set correctly.", i);
			}
		  locPot->AddValue(x1, x2, y);
		}
	  locPot->RecomputeLegendreRules();
	  potential = locPot;
	}
  else if(strcmp(temp.c_str(), "Parametrized") == 0)
	{
	  int potPrec;
	  if( ! poten.lookupValue("Precision", potPrec))
		{
		  throw RLException("Could not find precision setting in settings file.");
		}
	  Setting & pparam = poten["Parameters"];

	  vector<pair<string, double> > parameters;
	  for(int i = 0; i<pparam.getLength(); ++i)
		{
		  string paraStr;
		  double paraDbl;
		  if(!pparam[i].lookupValue("Name", paraStr) || !pparam[i].lookupValue("Value", paraDbl))
			{
			  throw RLException("Potential parameter #%d was not properly specified.", i);
			}
		  parameters.push_back(make_pair(paraStr, paraDbl));
		}
	  
	  string paraFunction;
	  if(!poten.lookupValue("Function", paraFunction))
		{
		  throw RLException("Could not find the 'Function' property in the parametrized potential in the config file.");
		}
	  double minX, maxX;
	  if(!poten.exists("Interval") || ! poten["Interval"].isArray() || poten["Interval"].getLength() != 2)
		{
		  throw RLException("Invalid interval specified for parametrized input potential.");
		}
	  minX = poten["Interval"][0];
	  maxX = poten["Interval"][1];
	  if(minX > maxX)
		{
		  throw RLException("Invalid interval specified in settings file: min is larger than max.");
		}
	  ParametrizedPotential * locPot = new ParametrizedPotential(paraFunction,
																 parameters, 
																 minX,
																 maxX
																 );
	  locPot->SetPrecision(potPrec);
	  potential = locPot;
	}
  else
	{
	  throw RLException("Unsupported potential type: '%s'", temp.c_str());
	}





  
  if(! computation.exists("KCurve") || !computation["KCurve"].isGroup())
	{
	  throw RLException("KCurve not properly set.");
	}
  
  Setting & kcur = computation["KCurve"];
  
  if( ! kcur.exists("Values") || ! kcur["Values"].isList())
	{
	  throw RLException("Values in KCurve not properly set.");
	}
  Setting & kv = kcur["Values"];
  ComplexDouble lastVal;
  for(int i = 0; i<kv.getLength(); ++i)
	{
	  if(!kv[i].exists("P0") || !kv[i]["P0"].isArray() || kv[i]["P0"].getLength() != 2)
		{
		  throw RLException("P0 not properly set in item %d.", i);
		}
	  if(!kv[i].exists("P1") || !kv[i]["P1"].isArray() || kv[i]["P1"].getLength() != 2)
		{
		  throw RLException("P1 not properly set in item %d.", i);
		}
	  double p0r, p0i, p1r, p1i;
	  p0r = kv[i]["P0"][0];
	  p0i = kv[i]["P0"][1];
	  p1r = kv[i]["P1"][0];
	  p1i = kv[i]["P1"][1];
	  ComplexDouble p0(p0r, p0i);
	  ComplexDouble p1(p1r, p1i);

	  int points;
	  if(! kv[i].lookupValue("Points", points) )
		{
		  throw RLException("Points not found in item %d.", i);
		}

	  if(i>0)
		{
		  if(!DBL_EQUAL(lastVal, p0))
			{
			  throw RLException("The specified curve had a discontinuity: %lf + i %lf is not equal to %lf + i %lf.",
								real(lastVal), imag(lastVal), p0r, p0i);
			}
		}
	  else
		{
		  kCurve->AddValue(p0);
		}
	  kCurve->AddGLPoints(points);
	  kCurve->AddValue(p1);
	  lastVal = p1;
	}
  
  kCurve->ComputeGaussLegendre();
  
  // Read the file. If there is an error, report it and exit.
}

unsigned int ComputeConfig::GetVerbosityLevel() const
{
  return verbosityLevel;
}

void ComputeConfig::SetVerbosityLevel(unsigned int value) 
{
  verbosityLevel = value;
}

void ComputeConfig::WriteFile(const char * fileName) const
{
  Config cfg;
  Setting &root = cfg.getRoot();
  
  root.add("Version", Setting::TypeFloat) = CONFIG_FILE_VERSION;
  root.add("Program", Setting::TypeGroup);

  root["Program"].add("VerbosityLevel", Setting::TypeInt) = (int)verbosityLevel;
  root["Program"].add("NumberOfThreads", Setting::TypeInt) = (int)numberOfThreads;

  Setting & autolaunch = root["Program"].add("AutoLaunch",Setting::TypeGroup);
  Setting & gnuplot = autolaunch.add("Gnuplot", Setting::TypeGroup);  
  gnuplot.add("Potential", Setting::TypeBoolean) = autoPlotPotential;
  gnuplot.add("KCurve", Setting::TypeBoolean) = autoPlotKCurve;

  Setting & output = root.add("Output", Setting::TypeGroup);
  output.add("KCurve", Setting::TypeString) = kCurveFile;
  output.add("KFound", Setting::TypeString) = kFoundFile;
  output.add("Potential", Setting::TypeString) = potentialFile;
  output.add("PotentialPrecision", Setting::TypeString) = potentialPrecisionFile;
  output.add("InterestingPoints", Setting::TypeString) = interestingPointsFile;
  output.add("Wavefunctions", Setting::TypeString) = wavefunctionFile;
  Setting & wavePrecision = output.add("WavefunctionProperties", Setting::TypeGroup);
  wavePrecision.add("MinX", Setting::TypeFloat) = minWavefunctionX;
  wavePrecision.add("MaxX", Setting::TypeFloat) = minWavefunctionX;
  wavePrecision.add("DeltaX", Setting::TypeFloat) = minWavefunctionX;
  
  
  root.add("Computation", Setting::TypeGroup);

  string matType = "General";
  if(matrixType == HermitianMatrix)
	matType = "Hermitian";
  if(matrixType == SymmetricMatrix)
	matType = "Symmetric";
  root["Computation"].add("ExpectedMatrixType", Setting::TypeString) = matType;
	

  Setting & basFun = root["Computation"].add("BasisFunctions", Setting::TypeArray);
  for(vector<BasisFunction>::const_iterator it = basisFunctions.begin(); it!=basisFunctions.end(); ++it)
	{
	  basFun.add(Setting::TypeString) = it->GetName();
	}

  Setting & poten = root["Computation"].add("Potential", Setting::TypeGroup);
  if(PiecewiseConstantPotential * locPot = dynamic_cast<PiecewiseConstantPotential*>(potential))
	{
	  poten.add("Type", Setting::TypeString) = "PiecewiseConstant";
	  poten.add("Precision", Setting::TypeInt) = 100;
	  Setting & values = poten.add("Values", Setting::TypeList);
	  
	  list<Interval> points = locPot->GetPotentialPoints();
	  for(list<Interval>::const_iterator it = points.begin(); it!= points.end(); ++it)
		{
		  Setting & p0 = values.add(Setting::TypeGroup);
		  Setting & r0 = p0.add("Interval", Setting::TypeArray);
		  r0.add(Setting::TypeFloat) = it->x1;
		  r0.add(Setting::TypeFloat) = it->x2;
		  p0.add("Value",Setting::TypeFloat) = it->y;
		}
	}

  Setting & kcurve = root["Computation"].add("KCurve", Setting::TypeGroup);
  Setting & kval = kcurve.add("Values", Setting::TypeList);


  for(unsigned int i = 0; i<kCurve->GetNumberOfSegments(); ++i)
	{
	  Setting & k0 = kval.add(Setting::TypeGroup);
	  Setting & k0p0 = k0.add("P0", Setting::TypeArray);
	  k0p0.add(Setting::TypeFloat) = real(kCurve->SegmentEvaluate(i,-1));
	  k0p0.add(Setting::TypeFloat) = imag(kCurve->SegmentEvaluate(i,-1));
	  Setting & k0p1 = k0.add("P1", Setting::TypeArray);
	  k0p1.add(Setting::TypeFloat) = real(kCurve->SegmentEvaluate(i,1));
	  k0p1.add(Setting::TypeFloat) = imag(kCurve->SegmentEvaluate(i,1));
	  k0.add("Points", Setting::TypeInt) = 20;
	}

  // Write out the updated configuration.
  try
  {
    cfg.writeFile(fileName);
  }
  catch(const FileIOException &fioex)
  {
	throw RLException("IO error when writing the file '%s'", fileName);
  }
  return;
}


double ComputeConfig::GetVersion() const
{
  return CONFIG_FILE_VERSION;
}

bool ComputeConfig::GetAutoPlotPotential() const
{
  return autoPlotPotential;
}

void ComputeConfig::SetAutoPlotPotential(bool value)
{
  autoPlotPotential = value;
}

bool ComputeConfig::GetAutoPlotKCurve() const
{
  return autoPlotKCurve;
}

void ComputeConfig::SetAutoPlotKCurve(bool value)
{
  autoPlotKCurve = value;
}

const char * ComputeConfig::GetKCurveFile() const
{
  return kCurveFile;
}

void ComputeConfig::SetKCurveFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(kCurveFile, value);
}



const char * ComputeConfig::GetKFoundFile() const
{
  return kFoundFile;
}

void ComputeConfig::SetKFoundFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(kFoundFile, value);
}


const char * ComputeConfig::GetPotentialFile() const
{
  return potentialFile;
}

void ComputeConfig::SetPotentialFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(potentialFile, value);
}

const char * ComputeConfig::GetPotentialPrecisionFile() const
{
  return potentialPrecisionFile;
}

void ComputeConfig::SetPotentialPrecisionFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(potentialPrecisionFile, value);
}

const char * ComputeConfig::GetInterestingPointsFile() const
{
  return interestingPointsFile;
}

void ComputeConfig::SetInterestingPointsFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(interestingPointsFile, value);
}

const char * ComputeConfig::GetWavefunctionFile() const
{
  return wavefunctionFile;
}

void ComputeConfig::SetWavefunctionFile(const char * value)
{
  if(strlen(value) > MAX_FILENAME_SIZE)
	throw RLException("Too long filename.");
  strcpy(wavefunctionFile, value);
}



Potential * ComputeConfig::GetPotential() const
{
  return potential;
}

void ComputeConfig::SetPotential(Potential * value)
{
  ///free stuff.
  if(potential != NULL)
	delete potential;

  potential = value;
}

ParametrizedCurve * ComputeConfig::GetKCurve() const
{
  return kCurve;
}

void ComputeConfig::SetKCurve(ParametrizedCurve * value)
{
  if(kCurve != NULL)
	delete kCurve;
  kCurve = value;
}

const vector<BasisFunction> & ComputeConfig::GetBasisFunctions() const
{
  return basisFunctions;
}

void ComputeConfig::SetBasisFunctions(vector<BasisFunction> value)
{
  basisFunctions = value;
}

ExpectedMatrixType ComputeConfig::GetExpectedMatrixType() const
{
  return matrixType;
}


void ComputeConfig::SetExpectedMatrixType(ExpectedMatrixType value)
{
  matrixType = value;
}


unsigned int ComputeConfig::GetNumberOfThreads() const
{
  return numberOfThreads;
}

void ComputeConfig::SetNumberOfThreads(unsigned int value)
{
  numberOfThreads = value;
}


void ComputeConfig::SetMinWavefunctionX(double value)
{
  minWavefunctionX = value;
}

void ComputeConfig::SetMaxWavefunctionX(double value)
{
  maxWavefunctionX = value;
}

void ComputeConfig::SetWavefunctionStepsizeX(double value)
{
  wavefunctionStepsizeX = value;
}


double ComputeConfig::GetMinWavefunctionX() const
{
  return minWavefunctionX;
}

double ComputeConfig::GetMaxWavefunctionX() const
{
  return maxWavefunctionX;
}

double ComputeConfig::GetWavefunctionStepsizeX() const
{
  return wavefunctionStepsizeX;
}
