#include "ComputeConfig.hh"

ComputeConfig::ComputeConfig()
  :verbosityLevel(0), autoPlotPotential(false), autoPlotKCurve(false)
{
  strcpy(kCurveFile, "KCurve.dat");
  strcpy(kFoundFile, "KFound.dat");
  strcpy(potentialFile,"Potential.dat");
  potential = new Potential();
  potential->AddValue(-3, -1, 10);
  potential->AddValue(-1, 1, -225);
  potential->AddValue(1, 3, 10);

  kCurve = new ParametrizedCurve(-1,1);
  kCurve->AddValue(0.0);
  kCurve->AddValue(ComplexDouble(3., 0.2));
  kCurve->AddValue(ComplexDouble(8.,0.));
  kCurve->AddValue(ComplexDouble(30.,0.));
  kCurve->AddGLPoints(25);
  kCurve->AddGLPoints(20);
  kCurve->AddGLPoints(30);
  kCurve->ComputeGaussLegendre();


  basisFunctions.push_back(BasisFunction("expi+"));
  basisFunctions.push_back(BasisFunction("expi-"));

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

  potential->Clear();
  kCurve->Clear();


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

  if( !root.exists("Computation") || !root["Computation"].isGroup())
	{
	  throw RLException("'Computation' was not appropriately defined as a group.");
	}

  Setting & computation = root["Computation"];
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
  if(strcmp(temp.c_str(), "PiecewiseConstant") != 0)
	{
	  throw RLException("Unsupported potential type: %s\n", temp.c_str());
	}
  
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
	  potential->AddValue(x1, x2, y);
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



void ComputeConfig::WriteFile(const char * fileName) const
{
  Config cfg;
  Setting &root = cfg.getRoot();
  
  root.add("Version", Setting::TypeFloat) = CONFIG_FILE_VERSION;
  root.add("Program", Setting::TypeGroup);

  root["Program"].add("VerbosityLevel", Setting::TypeInt) = (int)verbosityLevel;

  Setting & autolaunch = root["Program"].add("AutoLaunch",Setting::TypeGroup);
  Setting & gnuplot = autolaunch.add("Gnuplot", Setting::TypeGroup);  
  gnuplot.add("Potential", Setting::TypeBoolean) = autoPlotPotential;
  gnuplot.add("KCurve", Setting::TypeBoolean) = autoPlotKCurve;

  Setting & output = root.add("Output", Setting::TypeGroup);
  output.add("KCurve", Setting::TypeString) = kCurveFile;
  output.add("KFound", Setting::TypeString) = kFoundFile;
  output.add("Potential", Setting::TypeString) = potentialFile;
  
  root.add("Computation", Setting::TypeGroup);

  Setting & basFun = root["Computation"].add("BasisFunctions", Setting::TypeArray);
  for(list<BasisFunction>::const_iterator it = basisFunctions.begin(); it!=basisFunctions.end(); ++it)
	{
	  basFun.add(Setting::TypeString) = it->GetName();
	}

  Setting & poten = root["Computation"].add("Potential", Setting::TypeGroup);
  poten.add("Type", Setting::TypeString) = "PiecewiseConstant";
  Setting & values = poten.add("Values", Setting::TypeList);

  list<Interval> points = potential->GetPotentialPoints();
  for(list<Interval>::const_iterator it = points.begin(); it!= points.end(); ++it)
	{
	  Setting & p0 = values.add(Setting::TypeGroup);
	  Setting & r0 = p0.add("Interval", Setting::TypeArray);
	  r0.add(Setting::TypeFloat) = it->x1;
	  r0.add(Setting::TypeFloat) = it->x2;
	  p0.add("Value",Setting::TypeFloat) = it->y;
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
