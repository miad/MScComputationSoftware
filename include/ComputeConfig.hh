#ifndef ComputeConfig_hh
#define ComputeConfig_hh 1

#include "libconfig.hh"
#include "RLException.hh"

using namespace std;
using namespace libconfig;

#define CONFIG_FILE_VERSION 0.1

class ComputeConfig
{
public:
  ComputeConfig(); ///Constructor.
  ~ComputeConfig(); ///Destructor.
  void ReadFile(const char * fileName ///File name to read.
				); ///Read configuration from the specified file. Throws an exception if there is any problem with the parsing.
  void WriteFile(const char * fileName ///File name to write to.
				 ) const; ///Writes the settings to the specified file. Throws an exception if there is any problem.
private:

};



#endif
