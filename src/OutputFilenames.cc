#include "OutputFilenames.hh"

OutputFilenames::OutputFilenames()
{

}

void OutputFilenames::Add(string key, string value)
{
  filenames[key] = value;
}

void OutputFilenames::Clear()
{
  filenames.clear();
}

string OutputFilenames::Get(string key) const
{
  map<string, string>::const_iterator it = filenames.find(key);
  if(it == filenames.end())
	throw RLException("Attempted to access output filename for key '%s', which was not possible.", key.c_str());
  return it->second;
}
