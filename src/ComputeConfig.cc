#include "ComputeConfig.hh"

ComputeConfig::ComputeConfig()
  :numberOfThreads(1), verbosityLevel(0), autoPlotPotential(false), autoPlotKCurve(false),
   autoPlotWavefunctions(false),
   minWavefunctionX(-10), maxWavefunctionX(10), wavefunctionStepsizeX(0.01),
   numberOfParticles(2), couplingCoefficient(0.0),
   harmonicOverride(false), harmonicAngularFrequency(0.0), harmonicNmax(1)
{
  outputFilenames.Add("KCurveFile", "KCurve.dat");
  outputFilenames.Add("KFoundFile", "KFound.dat");
  outputFilenames.Add("PotentialFile", "Potential.dat");
  outputFilenames.Add("PotentialPrecisionFile", "PotentialPrecision.dat");
  outputFilenames.Add("InterestingPointsFile", "InterestingPoints.dat");
  outputFilenames.Add("WavefunctionFile", "Wavefunctions.dat");
  outputFilenames.Add("Matrix", "Matrix.dat");


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
}

void ComputeConfig::ReadFile(const char * fileName)
{
  Config cfg;
  try
	{
	  cfg.readFile(fileName);
	}
  catch(ParseException &ex)
	{
	  printf("Caught ParseException : %s Location : line %d in the config file.\n", ex.getError(), ex.getLine());
	  throw ex;
	}
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


  if(! program.lookupValue("VerbosityLevel", verbosityLevel))
	{
	  throw RLException("The value 'VerbosityLevel' was not defined.");
	}

  if(! program.lookupValue("NumberOfThreads", numberOfThreads))
	{
	  throw RLException("The value 'NumberOfThreads' was not defined.");
	}


  ReadAutoLaunch(program);


  if( !root.exists("Output") || ! root["Output"].isGroup())
	{
	  throw RLException("The value 'Output' was not set appropriately.");
	}

  Setting & output = root["Output"];

  ReadOutputFiles(output);

  if( !root.exists("OutputSpecifics") || ! root["OutputSpecifics"].isGroup())
	{
	  throw RLException("Could not find OutputSpecifics group.");
	}

  Setting & outputSpecifics = root["OutputSpecifics"];


  ReadWavefunctionProperties(outputSpecifics);
  ReadExtraInteresting(outputSpecifics);






  if( !root.exists("Computation") || !root["Computation"].isGroup())
	{
	  throw RLException("'Computation' was not appropriately defined as a group.");
	}

  Setting & computation = root["Computation"];


  ReadSpecificUnits(computation);
  ReadHarmonicOscillator(computation);
  ReadExpectedMatrixType(computation);

  ReadBasisFunctions(computation);
  ReadMultiParticleData(computation);
  ReadPotential(computation);
  ReadKCurve(computation);
  // Read the file. If there is an error, report it and exit.
}

void ComputeConfig::ReadMultiParticleData(Setting & computation)
{
  if (! computation.exists("NumberOfParticles"))
	{
	  throw RLException("NumberOfParticles not properly specified.");
	}
  if(! computation.lookupValue("NumberOfParticles", numberOfParticles))
	{
	  throw RLException("Could not lookup # particles.");
	}
  if(numberOfParticles < 1 || numberOfParticles > 2)
	{
	  throw RLException("Unsupported # particles.");
	}

  if (! computation.exists("CouplingCoefficient"))
	{
	  throw RLException("CouplingCoefficient not properly specified.");
	}
  if(! computation.lookupValue("CouplingCoefficient", couplingCoefficient))
	{
	  throw RLException("Could not lookup coupling coefficient.");
	}

}

void ComputeConfig::ReadWavefunctionProperties(Setting & outputSpecifics)
{
  if( !outputSpecifics.exists("WavefunctionProperties") || !outputSpecifics["WavefunctionProperties"].isGroup() )
	{
	  throw RLException("WavefunctionProperties not properly defined.");
	}

  Setting & wfProp = outputSpecifics["WavefunctionProperties"];
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
}

void ComputeConfig::ReadExtraInteresting(Setting & outputSpecifics)
{
  extraInterestingPoints.clear();
  if( ! outputSpecifics.exists("ExtraInteresting") || ! outputSpecifics["ExtraInteresting"].isList() )
	{
	  throw RLException("ExtraInteresting not properly specified.");
	}
  Setting & extraInteresting = outputSpecifics["ExtraInteresting"];
  for(int i = 0; i<extraInteresting.getLength(); ++i)
	{
	  if(!extraInteresting[i].isArray() || (extraInteresting[i].getLength() != 2))
		{
		  throw RLException("Invalid number specification at number %d in ExtraInteresting list.", i);
		}
	  for(int j = 0; j<2; ++j)
		if(extraInteresting[i][j].getType() != Setting::TypeFloat )
		  {
			throw RLException("ExtraInteresting number %d position %d was not float.", i, j);
		  }
	  extraInterestingPoints.push_back(ComplexDouble(extraInteresting[i][0], extraInteresting[i][1]));
	}

  
}


uint ComputeConfig::GetVerbosityLevel() const
{
  return verbosityLevel;
}

void ComputeConfig::SetVerbosityLevel(uint value) 
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
  gnuplot.add("Wavefunctions", Setting::TypeBoolean) = autoPlotWavefunctions;

  Setting & output = root.add("Output", Setting::TypeGroup);
  output.add("KCurve", Setting::TypeString) = outputFilenames.Get("KCurveFile");
  output.add("KFound", Setting::TypeString) = outputFilenames.Get("KFoundFile");
  output.add("Potential", Setting::TypeString) = outputFilenames.Get("PotentialFile");
  output.add("PotentialPrecision", Setting::TypeString) = outputFilenames.Get("PotentialPrecisionFile");
  output.add("InterestingPoints", Setting::TypeString) = outputFilenames.Get("InterestingPointsFile");
  output.add("Wavefunctions", Setting::TypeString) = outputFilenames.Get("WavefunctionFile");
  output.add("Matrix", Setting::TypeString) = outputFilenames.Get("MatrixFile");
  Setting & outputSpecifics = root.add("OutputSpecifics", Setting::TypeGroup);
  Setting & wavePrecision = outputSpecifics.add("WavefunctionProperties", Setting::TypeGroup);
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
  root["Computation"].add("NumberOfParticles", Setting::TypeInt) = (int)numberOfParticles;
  root["Computation"].add("CouplingCoefficient", Setting::TypeFloat) = couplingCoefficient;
	

  Setting & units = root["Computation"].add("Units", Setting::TypeGroup);
  units.add("HbarTimesLambda", Setting::TypeFloat) = specificUnits.GetHbarTimesLambda();
  units.add("MassOverLambda2", Setting::TypeFloat) = specificUnits.GetMassOverLambda2();
  units.add("LengthUnitName", Setting::TypeString) = specificUnits.GetLengthUnitName();
  units.add("EnergyUnitName", Setting::TypeString) = specificUnits.GetEnergyUnitName();

  Setting & oscillator = root["Computation"].add("HarmonicOscillator", Setting::TypeGroup);
  oscillator.add("Override", Setting::TypeBoolean) = harmonicOverride;
  oscillator.add("AngularFrequency", Setting::TypeFloat) = harmonicAngularFrequency;
  oscillator.add("Nmax", Setting::TypeInt) = (int)harmonicNmax;


  Setting & basFun = root["Computation"].add("BasisFunctions", Setting::TypeArray);
  for(vector<BasisFunction>::const_iterator it = basisFunctions.begin(); it!=basisFunctions.end(); ++it)
	{
	  basFun.add(Setting::TypeString) = it->GetName();
	}

  Setting & poten = root["Computation"].add("PotentialFile", Setting::TypeGroup);
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


  for(uint i = 0; i<kCurve->GetNumberOfSegments(); ++i)
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

void ComputeConfig::SetExtraInterestingPoints(vector<ComplexDouble> value)
{
  extraInterestingPoints = value;
}

const vector<ComplexDouble> & ComputeConfig::GetExtraInterestingPoints() const
{
  return extraInterestingPoints;
}

bool ComputeConfig::GetAutoPlotKCurve() const
{
  return autoPlotKCurve;
}

void ComputeConfig::SetAutoPlotKCurve(bool value)
{
  autoPlotKCurve = value;
}

bool ComputeConfig::GetAutoPlotWavefunctions() const
{
  return autoPlotWavefunctions;
}

void ComputeConfig::SetAutoPlotWavefunctions(bool value)
{
  autoPlotWavefunctions = value;
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


uint ComputeConfig::GetNumberOfThreads() const
{
  return numberOfThreads;
}

void ComputeConfig::SetNumberOfThreads(uint value)
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


const OutputFilenames * ComputeConfig::GetOutputFilenames() const
{
  return &outputFilenames;
}

void ComputeConfig::SetOutputFilenames(const OutputFilenames & value)
{
  outputFilenames = value;
}

const SpecificUnits * ComputeConfig::GetSpecificUnits() const
{
  return &specificUnits;
}

void ComputeConfig::SetSpecificUnits(const SpecificUnits & value)
{
  specificUnits = value;
}

void ComputeConfig::ReadOutputFiles(Setting & output)
{
  outputFilenames.Clear();
  string temp;

  if( output.lookupValue("KCurve", temp) )
	{
	  outputFilenames.Add("KCurveFile", temp);
	}
  else
	{
	  throw RLException("The value 'KCurve' was not set appropriately.");
	}


  if( output.lookupValue("KFound", temp) )
	{
	  outputFilenames.Add("KFoundFile", temp);
	}
  else
	{
	  throw RLException("The value 'KFound' was not set appropriately.");
	}


  if( output.lookupValue("Potential", temp) )
	{
	  outputFilenames.Add("PotentialFile", temp);
	}
  else
	{
	  throw RLException("The value 'Potential' was not set appropriately.");
	}


  if( output.lookupValue("PotentialPrecision", temp) )
	{
	  outputFilenames.Add("PotentialPrecisionFile", temp);
	}
  else
	{
	  throw RLException("The value 'PotentialPrecision' was not set appropriately.");
	}


  if( output.lookupValue("InterestingPoints", temp) )
	{
	  outputFilenames.Add("InterestingPointsFile", temp);
	}
  else
	{
	  throw RLException("The value 'InterestingPoints' was not set appropriately.");
	}


  if( output.lookupValue("Wavefunctions", temp) )
	{
	  outputFilenames.Add("WavefunctionsFile", temp);
	}
  else
	{
	  throw RLException("The value 'Wavefunctions' was not set appropriately.");
	}

  if( output.lookupValue("Matrix", temp) )
	{
	  outputFilenames.Add("MatrixFile", temp);
	}
  else
	{
	  throw RLException("The value 'Matrix' was not set appropriately.");
	}

}


void ComputeConfig::ReadExpectedMatrixType(Setting & computation)
{
  string temp;

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
}

void ComputeConfig::ReadSpecificUnits(Setting & computation)
{
  if (! computation.exists("Units") || !computation["Units"].isGroup())
	{
	  throw RLException("Units group not properly specified in settings file.");
	}
  
  Setting & units = computation["Units"];
  double hbarTimesLambda, massOverLambda2;
  if( ! units.lookupValue("HbarTimesLambda", hbarTimesLambda) )
	{
	  throw RLException("HbarTimesLambda unit not properly specified.");
	}
  if( ! units.lookupValue("MassOverLambda2", massOverLambda2) )
	{
	  throw RLException("MassOverLambda2 unit not properly specified.");
	}
  string lengthUnitName, energyUnitName;
  if( ! units.lookupValue("LengthUnitName", lengthUnitName) )
	{
	  throw RLException("Could not find LengthUnitName in config file.");
	}
  if( !units.lookupValue("EnergyUnitName", energyUnitName) )
	{
	  throw RLException("Could not find EnergyUnitName in config file.");
	}
  specificUnits = SpecificUnits(hbarTimesLambda, massOverLambda2, lengthUnitName, energyUnitName);
}

void ComputeConfig::ReadHarmonicOscillator(Setting & computation)
{
  if (! computation.exists("HarmonicOscillator") || !computation["HarmonicOscillator"].isGroup())
	{
	  throw RLException("HarmonicOscillator group not properly specified in settings file.");
	}
  
  Setting & oscillator = computation["HarmonicOscillator"];

  if( ! oscillator.lookupValue("Override", harmonicOverride) )
	{
	  throw RLException("Could not locate harmonic override.");
	}
  if( ! oscillator.lookupValue("AngularFrequency", harmonicAngularFrequency ) )
	{
	  throw RLException("Could not locate harmonic angular frequency.");
	}
  if( ! oscillator.lookupValue("Nmax", harmonicNmax) )
	{
	  throw RLException("Could not locate harmonicNmax");
	}
}



void ComputeConfig::ReadPotential(Setting & computation)
{

  if( !computation.exists("Potential") || !computation["Potential"].isGroup())
	{
	  throw RLException("Potential not properly defined.");
	}

  Setting & poten = computation["Potential"];

  string temp;
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

	  if( ! poten.exists("Parameters") || ! poten["Parameters"].isList() )
		{
		  throw RLException("Potential parameters not properly specified.");
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
	  string minX, maxX;
	  if(!poten.exists("Interval") || ! poten["Interval"].isArray() || poten["Interval"].getLength() != 2)
		{
		  throw RLException("Invalid interval specified for parametrized input potential.");
		}
	  
	  minX = poten["Interval"][0].c_str();
	  maxX = poten["Interval"][1].c_str();

	  ParametrizedPotential * locPot = new ParametrizedPotential(paraFunction,
																 parameters, 
																 minX.c_str(),
																 maxX.c_str()
																 );
	  locPot->SetPrecision(potPrec);
	  potential = locPot;
	}
  else
	{
	  throw RLException("Unsupported potential type: '%s'", temp.c_str());
	}
}


void ComputeConfig::ReadKCurve(Setting & computation)
{

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
  

}


void ComputeConfig::ReadBasisFunctions(Setting & computation)
{
  if( !computation.exists("BasisFunctions") || !computation["BasisFunctions"].isArray())
	{
	  throw RLException("Basis functions not properly defined.");
	}

  Setting & bf = computation["BasisFunctions"];
  for(int i = 0; i<bf.getLength(); ++i)
	{
	  basisFunctions.push_back(BasisFunction(bf[i]));
	}
}


void ComputeConfig::ReadAutoLaunch(Setting & program)
{
  if(!program.exists("AutoLaunch") || ! program["AutoLaunch"].isGroup())
	{
	  throw RLException("The group 'AutoLaunch' was not properly defined in the config file.");
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

  if( ! gnuplot.lookupValue("Wavefunctions", autoPlotWavefunctions))
	{
	  throw RLException("The value 'Wavefunctions' in 'AutoLaunch' was not set appropriately.");
	}
}

void ComputeConfig::SetNumberOfParticles(uint value)
{
  numberOfParticles = value;
}

void ComputeConfig::SetCouplingCoefficient(double value)
{
  couplingCoefficient = value;
}

uint ComputeConfig::GetNumberOfParticles() const
{
  return numberOfParticles;
}

double ComputeConfig::GetCouplingCoefficient() const
{
  return couplingCoefficient;
}

bool ComputeConfig::GetHarmonicOverride() const
{
  return harmonicOverride;
}

double ComputeConfig::GetHarmonicAngularFrequency() const
{
  return harmonicAngularFrequency;
}

uint ComputeConfig::GetHarmonicNmax() const
{
  return harmonicNmax;
}

void ComputeConfig::SetHarmonicOverride(bool value)
{
  harmonicOverride = value;
}

void ComputeConfig::SetHarmonicAngularFrequency(double value)
{
  harmonicAngularFrequency = value;
}

void ComputeConfig::SetHarmonicNmax(uint value)
{
  harmonicNmax = value;
}
