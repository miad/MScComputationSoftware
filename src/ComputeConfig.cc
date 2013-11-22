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
  
  // Read the file. If there is an error, report it and exit.
}

void ComputeConfig::WriteFile(const char * fileName) const
{
  Config cfg;
  Setting &root = cfg.getRoot();
  
  if(! root.exists("version"))
    root.add("version", Setting::TypeFloat) = CONFIG_FILE_VERSION;

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
