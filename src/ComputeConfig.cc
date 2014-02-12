#include "ComputeConfig.hh"

ComputeConfig::ComputeConfig()
  :numberOfThreads(1), verbosityLevel(0), autoPlotPotential(false), autoPlotKCurve(false),
   autoPlotWavefunctions(false),
   minWavefunctionX(-10), maxWavefunctionX(10), wavefunctionStepsizeX(0.01),
   numberOfParticles(2),
   harmonicOverride(false),harmonicNmax(1)
{
  outputFilenames.Add("KCurveFile", "KCurve.dat");
  outputFilenames.Add("KFoundFile", "KFound.dat");
  outputFilenames.Add("PotentialFile", "Potential.dat");
  outputFilenames.Add("PotentialPrecisionFile", "PotentialPrecision.dat");
  outputFilenames.Add("InterestingPointsFile", "InterestingPoints.dat");
  outputFilenames.Add("WavefunctionFile", "Wavefunctions.dat");
  outputFilenames.Add("Matrix", "Matrix.dat");
  outputFilenames.Add("EnergyFile", "Energies.dat");
  outputFilenames.Add("ProductWavefunctionFile", "ProductWavefunction.dat");

  kCurve = new ParametrizedCurve(-1,1);
  kCurve->AddValue(0.0);
  kCurve->AddValue(ComplexDouble(3., -0.2));
  kCurve->AddValue(ComplexDouble(8.,0.));
  kCurve->AddValue(ComplexDouble(30.,0.));
  kCurve->AddGLPoints(25);
  kCurve->AddGLPoints(20);
  kCurve->AddGLPoints(30);
  kCurve->ComputeGaussLegendre();

  myHarmonicBasisFunction = NULL; ///Note: this will seriously impede functionality. Need to be fixed.
	//new HarmonicBasisFunction(vector<double>(2, 0), vector<double>(2, 1), &potentials, &specificUnits);

  basisFunctions.push_back(BasisFunction("sqrt(2)*sin(k*x)"));
  basisFunctions.push_back(BasisFunction("sqrt(2)*cos(k*x)"));
}

ComputeConfig::ComputeConfig(const char * fileName)
{
  ReadFile(fileName);
}

ComputeConfig::~ComputeConfig()
{
  for(vector<Potential*>::iterator it = potentials.begin(); it!=potentials.end(); ++it)
	{
	  if(*it != NULL)
		delete *it;
	}
  potentials.clear();


  if(kCurve != NULL)
	{
	  delete kCurve;
	  kCurve = NULL;
	}
  if(myHarmonicBasisFunction != NULL)
	{
	  delete myHarmonicBasisFunction;
	  myHarmonicBasisFunction = NULL;
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

  for(vector<Potential*>::iterator it = potentials.begin(); it!=potentials.end(); ++it)
	{
	  delete *it;
	}
  potentials.clear();


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
  ReadExpectedMatrixType(computation);

  ReadBasisFunctions(computation);
  ReadMultiParticleData(computation);
  ReadInteractionProperties(computation);
  ReadPotential(computation);
  ReadKCurve(computation);

  ReadHarmonicOscillator(computation);

  ReadSolverInfo(computation);

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
}

void ComputeConfig::ReadInteractionProperties(Setting & computation)
{
  if( ! computation.exists("Interaction") || ! computation["Interaction"].isGroup() )
	{
	  throw RLException("Interaction group not properly specified in config file.");
	}
  Setting & inter = computation["Interaction"];

  if (! inter.exists("CouplingCoefficient"))
	{
	  throw RLException("CouplingCoefficient not properly specified.");
	}
  double cpCoeff;
  if(! inter.lookupValue("CouplingCoefficient", cpCoeff))
	{
	  throw RLException("Could not lookup coupling coefficient.");
	}
  myInteractionProperties.SetCouplingCoefficient(cpCoeff);

  int nmax;
  if(! inter.lookupValue("Nmax", nmax) )
	{
	  throw RLException("Could not look up nmax for interaction.");
	}
  myInteractionProperties.SetNMax(nmax);

  int prec;
  if(! inter.lookupValue("Precision", prec) )
	{
	  throw RLException("Could not find interaction precision in config file.");
	}
  myInteractionProperties.SetPrecision(prec);

  string cFile;
  if(! inter.lookupValue("CacheFile", cFile) )
	{
	  throw RLException("Could not find cache file in config file.");
	}
  myInteractionProperties.SetCacheFile(cFile);

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
  output.add("EnergyFile", Setting::TypeString) = outputFilenames.Get("EnergyFile");
  output.add("ProductWavefunction", Setting::TypeString) = outputFilenames.Get("ProductWavefunction");
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
  Setting & inter = root["Computation"].add("Interaction", Setting::TypeGroup);
  inter.add("CouplingCoefficient", Setting::TypeFloat) = myInteractionProperties.GetCouplingCoefficient();
  inter.add("Nmax", Setting::TypeInt) = (int)myInteractionProperties.GetNMax(); 
  inter.add(Setting::TypeInt) = (int)myInteractionProperties.GetPrecision();
  inter.add("CacheFile", Setting::TypeString) = myInteractionProperties.GetCacheFile();


  Setting & locSolver = root["Computation"].add("Solver", Setting::TypeGroup);


  string sType = "Lapack";
  if(mySolver.GetSolverType() == ArpackSolver)
	sType = "Arpack";
  locSolver.add("Type", Setting::TypeString) = sType;
  locSolver.add("WorkArea", Setting::TypeString) = mySolver.GetWorkArea();
  locSolver.add("NumberOfValues", Setting::TypeInt) = (int)mySolver.GetNumberOfEigenvalues();
  Setting & shiftArray = locSolver.add("Shift", Setting::TypeArray);
  shiftArray.add( Setting::TypeFloat) = real(mySolver.GetShift());
  shiftArray.add( Setting::TypeFloat) = imag(mySolver.GetShift());
  
	

  Setting & units = root["Computation"].add("Units", Setting::TypeGroup);
  units.add("HbarTimesLambda", Setting::TypeFloat) = specificUnits.GetHbarTimesLambda();
  units.add("MassOverLambda2", Setting::TypeFloat) = specificUnits.GetMassOverLambda2();
  units.add("LengthUnitName", Setting::TypeString) = specificUnits.GetLengthUnitName();
  units.add("EnergyUnitName", Setting::TypeString) = specificUnits.GetEnergyUnitName();
  units.add("TimeToHertzFactor", Setting::TypeFloat) = specificUnits.GetTimeToHertzFactor();

  Setting & oscillator = root["Computation"].add("HarmonicOscillator", Setting::TypeGroup);
  oscillator.add("Override", Setting::TypeBoolean) = harmonicOverride;
  oscillator.add("Nmax", Setting::TypeInt) = (int)harmonicNmax;

  Setting & oscillatorArray = oscillator.add("AngularFrequency", Setting::TypeArray);
  Setting & oscillatorXmin = oscillator.add("Xmin", Setting::TypeArray);
  for(uint i = 0; i<2; ++i)
	{
	  oscillatorArray.add(Setting::TypeFloat) = myHarmonicBasisFunction->GetOmega(i);
	  oscillatorXmin.add(Setting::TypeFloat) = myHarmonicBasisFunction->GetXmin(i);
	}



  Setting & basFun = root["Computation"].add("BasisFunctions", Setting::TypeArray);
  for(vector<BasisFunction>::const_iterator it = basisFunctions.begin(); it!=basisFunctions.end(); ++it)
	{
	  basFun.add(Setting::TypeString) = it->GetName();
	}

  Setting & poten = root["Computation"].add("PotentialFile", Setting::TypeGroup);
  if( potentials.size() > 0 )
	{
	  if(PiecewiseConstantPotential * locPot = dynamic_cast<PiecewiseConstantPotential*>(potentials.at(0)))
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


Potential * ComputeConfig::GetPotential(uint index) const
{
  if(index >= potentials.size())
	throw RLException("Invalid potential index: %d\n", index);
  return potentials[index];
}

void ComputeConfig::SetPotential(Potential * value, uint index)
{
  if(index < potentials.size())
	{
	  if(potentials.at(index) != NULL)
		{
		  delete potentials.at(index);
		}
	  potentials.at(index) = value;
	}
  else if(index == potentials.size() && numberOfParticles < potentials.size())
	{
	  potentials.push_back(value);
	}
  else
	{
	  throw RLException("Tried to set invalid potential index: %d\n", index);
	}
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

  if( output.lookupValue("EnergyFile", temp) )
	{
	  outputFilenames.Add("EnergyFile", temp);
	}
  else
	{
	  throw RLException("The value 'EnergyFile' was not set appropriately.");
	}

  if( output.lookupValue("ProductWavefunction", temp) )
	{
	  outputFilenames.Add("ProductWavefunction", temp);
	}
   else
	{
	  throw RLException("The value 'ProductWavefunction' was not set appropriately.");
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
  double hbarTimesLambda, massOverLambda2, timeToHertzFactor;
  if( ! units.lookupValue("HbarTimesLambda", hbarTimesLambda) )
	{
	  throw RLException("HbarTimesLambda unit not properly specified.");
	}
  if( ! units.lookupValue("MassOverLambda2", massOverLambda2) )
	{
	  throw RLException("MassOverLambda2 unit not properly specified.");
	}
  if(! units.lookupValue("TimeToHertzFactor", timeToHertzFactor) )
	{
	  throw RLException("TimeToHertzFactor undefined.");
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
  specificUnits = SpecificUnits(hbarTimesLambda, massOverLambda2, lengthUnitName, energyUnitName, timeToHertzFactor);
}

void ComputeConfig::ReadSolverInfo(Setting & computation)
{

  if( ! computation.exists("Solver") || ! computation["Solver"].isGroup())
	{
	  throw RLException("No info on eigenvalue solver found in config file.");
	}

  Setting & solver = computation["Solver"];
  string temp;
  if(! solver.lookupValue("Type", temp) )
	{
	  throw RLException("Could not look up value for Solver::Type.");
	}
  bool found = false;

  EigenSolverType mySolverType;
  if(strcmp(temp.c_str(), "Lapack") == 0) 
	{
	  mySolverType = LapackSolver;
	  found = true;
	}
  if(strcmp(temp.c_str(), "Arpack") == 0)
	{
	  mySolverType = ArpackSolver;
	  found = true;
	}
  mySolver.SetSolverType(mySolverType);

  if(!found)
	throw RLException("Invalid option for Solver::Type.");

  if(! solver.lookupValue("WorkArea", temp))
	{
	  throw RLException("Could not look up value for Solver::WorkArea.");
	}
  mySolver.SetWorkArea(temp);
  if(! solver.exists("Shift") || ! solver["Shift"].isArray() || solver["Shift"].getLength() != 2)
	{
	  throw RLException("Invalid solver shift.");
	}
  mySolver.SetShift(ComplexDouble(solver["Shift"][0], solver["Shift"][1]));
  uint temp2;
  if(!solver.lookupValue("NumberOfValues", temp2))
	{
	  throw RLException("Could not look up NumberOfValues for solver.");
	}
  mySolver.SetNumberOfEigenvalues(temp2);
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
  if( ! oscillator.lookupValue("Nmax", harmonicNmax) )
	{
	  throw RLException("Could not locate harmonicNmax");
	}
  vector<double> harmonicAngularFrequency, harmonicXMin;
  
  if( ! oscillator.exists("AngularFrequency") || ! oscillator["AngularFrequency"].isArray() || oscillator["AngularFrequency"].getLength() != 2)
	{
	  throw RLException("Invalid array at AngularFrequency: should contain two floats.");
	}
  harmonicAngularFrequency.push_back(oscillator["AngularFrequency"][0]);
  harmonicAngularFrequency.push_back(oscillator["AngularFrequency"][1]);


  if( ! oscillator.exists("Xmin") || ! oscillator["Xmin"].isArray() || oscillator["Xmin"].getLength() != 2)
	{
	  throw RLException("Invalid array at XMin: should contain two floats.");
	}
  harmonicXMin.push_back(oscillator["Xmin"][0]);
  harmonicXMin.push_back(oscillator["Xmin"][1]);

  if(myHarmonicBasisFunction)
	{
	  delete myHarmonicBasisFunction;
	}

  myHarmonicBasisFunction = new HarmonicBasisFunction(harmonicXMin, harmonicAngularFrequency, &potentials, &specificUnits, potentials.at(0)->GetPrecision());
}

void ComputeConfig::ReadPiecewiseConstantPotential(Setting & poten)
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
  potentials.push_back(locPot);
}

void ComputeConfig::ReadParametrizedPotential(Setting & poten)
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



  map<string, vector<double> > variableParameters;

  ///Var-param will superseed corresponding non-var-params.
  if( ! poten.exists("VariableParameters") || ! poten["VariableParameters"].isList() )
	{
	  throw RLException("VariableParameters for potential not properly specified.");
	}
  Setting & vparam = poten["VariableParameters"];
  
  for(int i = 0; i<vparam.getLength(); ++i)
	{
	  string paraStr;
	  if(!vparam[i].lookupValue("Name", paraStr))
		{
		  throw RLException("Variable parameter name #%d was not properly specified.", i);
		}
	  variableParameters[paraStr] = vector<double>();
	  if(! vparam[i].exists("Value") || ! vparam[i]["Value"].isArray() || vparam[i]["Value"].getLength() != (int)2 )
		{
		  throw RLException("Invalid value given for variable parameter '%s'.", paraStr.c_str());
		}
	  for(int j = 0; j<vparam[i]["Value"].getLength(); ++j)
		{
		  variableParameters[paraStr].push_back(vparam[i]["Value"][j]);
		}
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
	  if(variableParameters.find(paraStr) == variableParameters.end())
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
  

  for(uint i = 0; i<2; ++i)
	{
	  vector<pair<string, double> > fullParam = parameters;
	  for(map<string, vector<double> >::const_iterator it = variableParameters.begin(); it!=variableParameters.end(); ++it)
		{
		  fullParam.push_back(make_pair(it->first, it->second[i]));
		}

	  ParametrizedPotential * locPot = new ParametrizedPotential(paraFunction,
																 fullParam, 
																 minX.c_str(),
																 maxX.c_str()
																 );
	  locPot->SetPrecision(potPrec);
	  potentials.push_back(locPot);
	}
}



void ComputeConfig::ReadPotential(Setting & computation)
{
  for(vector<Potential*>::iterator it = potentials.begin(); it!=potentials.end(); ++it)
	{
	  delete *it;
	}
  potentials.clear();



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
	  ReadPiecewiseConstantPotential(poten);
	}
  else if(strcmp(temp.c_str(), "Parametrized") == 0)
	{
	  ReadParametrizedPotential(poten);
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

void ComputeConfig::SetInteractionProperties(InteractionProperties value)
{
  myInteractionProperties = value;
}


uint ComputeConfig::GetNumberOfParticles() const
{
  return numberOfParticles;
}

const InteractionProperties * ComputeConfig::GetInteractionProperties() const
{
  return &myInteractionProperties;
}

bool ComputeConfig::GetHarmonicOverride() const
{
  return harmonicOverride;
}

uint ComputeConfig::GetHarmonicNmax() const
{
  return harmonicNmax;
}

void ComputeConfig::SetHarmonicOverride(bool value)
{
  harmonicOverride = value;
}

void ComputeConfig::SetHarmonicNmax(uint value)
{
  harmonicNmax = value;
}

HarmonicBasisFunction * ComputeConfig::GetHarmonicBasisFunction() const
{
  return myHarmonicBasisFunction;
}

void ComputeConfig::SetHarmonicBasisFunction(HarmonicBasisFunction * value)
{
  myHarmonicBasisFunction = value;
}

EigenvalueSolver * ComputeConfig::GetSolver()
{
  return &mySolver;
}
