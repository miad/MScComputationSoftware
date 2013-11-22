#include "ComputeConfig.hh"

ComputeConfig::ComputeConfig()
{

}

ComputeConfig::~ComputeConfig()
{

}

void ComputeConfig::ReadFile(const char * fileName)
{
  Config cfg;
  cfg.readFile(fileName);
  Setting & root = cfg.getRoot();


  ///Validate version.
  double version;
  if( root.lookupValue("Version", version) )
	{
	  if( ! DBL_EQUAL(version, CONFIG_FILE_VERSION) )
		throw RLException("Invalid configuration file version: was: %f, required %f", version, CONFIG_FILE_VERSION);
	}
  else
	{
	  throw RLException("Could not find property 'Version' in configuration file.");
	}

  


  // Read the file. If there is an error, report it and exit.
}

void ComputeConfig::WriteFile(const char * fileName) const
{
  Config cfg;
  Setting &root = cfg.getRoot();
  
  root.add("Version", Setting::TypeFloat) = CONFIG_FILE_VERSION;
  root.add("Program", Setting::TypeGroup);

  root["Program"].add("VerbosityLevel", Setting::TypeInt) = 0;

  Setting & autolaunch = root["Program"].add("AutoLaunch",Setting::TypeGroup);
  Setting & gnuplot = autolaunch.add("Gnuplot", Setting::TypeGroup);  
  gnuplot.add("Potential", Setting::TypeBoolean) = false;
  gnuplot.add("KCurve", Setting::TypeBoolean) = true;

  Setting & output = root.add("Output", Setting::TypeGroup);
  output.add("KCurve", Setting::TypeString) = "OutputK.dat";
  output.add("KFound", Setting::TypeString) = "OutputF.dat";
  output.add("Potential", Setting::TypeString) = "OutputV.dat";
  
  

  root.add("Computation", Setting::TypeGroup);

  Setting & potential = root["Computation"].add("Potential", Setting::TypeGroup);
  potential.add("Type", Setting::TypeString) = "PiecewiseConstant";
  Setting & values = potential.add("Values", Setting::TypeList);

  Setting & p0 = values.add(Setting::TypeGroup);
  Setting & r0 = p0.add("Interval", Setting::TypeArray);
  r0.add(Setting::TypeFloat) = -3.0;
  r0.add(Setting::TypeFloat) = -1.0;
  p0.add("Value",Setting::TypeFloat) = 10.0;

  Setting & p1 = values.add(Setting::TypeGroup);
  Setting & r1 = p1.add("Interval", Setting::TypeArray);
  r1.add(Setting::TypeFloat) = -1.0;
  r1.add(Setting::TypeFloat) = 1.0;
  p1.add("Value",Setting::TypeFloat) = -225.0;

  Setting & p2 = values.add(Setting::TypeGroup);
  Setting & r2 = p2.add("Interval", Setting::TypeArray);
  r2.add(Setting::TypeFloat) = 1.0;
  r2.add(Setting::TypeFloat) = 3.0;
  p2.add("Value",Setting::TypeFloat) = 10.0;


  Setting & kcurve = root["Computation"].add("KCurve", Setting::TypeGroup);
  Setting & kval = kcurve.add("Values", Setting::TypeList);
  
  Setting & k0 = kval.add(Setting::TypeGroup);
  Setting & k0p0 = k0.add("P0", Setting::TypeArray);
  k0p0.add(Setting::TypeFloat) = 0.0;
  k0p0.add(Setting::TypeFloat) = 0.0;
  Setting & k0p1 = k0.add("P1", Setting::TypeArray);
  k0p1.add(Setting::TypeFloat) = 3.0;
  k0p1.add(Setting::TypeFloat) = 0.2;
  k0.add("Points", Setting::TypeInt) = 20;

  Setting & k2 = kval.add(Setting::TypeGroup);
  Setting & k2p0 = k2.add("P0", Setting::TypeArray);
  k2p0.add(Setting::TypeFloat) = 3.0;
  k2p0.add(Setting::TypeFloat) = 0.2;
  Setting & k2p1 = k2.add("P1", Setting::TypeArray);
  k2p1.add(Setting::TypeFloat) = 8.0;
  k2p1.add(Setting::TypeFloat) = 0.0;
  k2.add("Points", Setting::TypeInt) = 20;


  Setting & k1 = kval.add(Setting::TypeGroup);
  Setting & k1p0 = k1.add("P0", Setting::TypeArray);
  k1p0.add(Setting::TypeFloat) = 8.0;
  k1p0.add(Setting::TypeFloat) = 0.0;
  Setting & k1p1 = k1.add("P1", Setting::TypeArray);
  k1p1.add(Setting::TypeFloat) = 30.0;
  k1p1.add(Setting::TypeFloat) = 0.0;
  k1.add("Points", Setting::TypeInt) = 30;
  

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
