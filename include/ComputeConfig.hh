#ifndef ComputeConfig_hh
#define ComputeConfig_hh 1

#include "libconfig.hh"
#include "RLException.hh"
#include "Globals.hpp"
#include "Potential.hh"
#include "ParametrizedCurve.hh"
#include "BasisFunction.hh"
#include <list>
#include "RLMacros.hpp"

using namespace std;
using namespace libconfig;

#define CONFIG_FILE_VERSION 0.1
#define MAX_FILENAME_SIZE 1000


  enum ExpectedMatrixType
	{
	  GeneralMatrix = 0,
	  SymmetricMatrix = 1,
	  HermitianMatrix = 2
	};




/* Configuration for the program.
   All options are \emph{guaranteed} to be valid.
 */
class ComputeConfig
{
public:
  ComputeConfig(); ///Constructor.
  ComputeConfig(const char * fileName ///Configuration file name.
				); ///Initialize directly from file. Faster than first constructing and then calling ReadFile.
  ~ComputeConfig(); ///Destructor.


  void ReadFile(const char * fileName ///File name to read.
				); ///Read configuration from the specified file. Throws an exception if there is any problem with the parsing.
  void WriteFile(const char * fileName ///File name to write to.
				 ) const; ///Writes the settings to the specified file. Throws an exception if there is any problem.



  ///Getters
  double GetVersion() const;

  unsigned int GetVerbosityLevel() const;

  bool GetAutoPlotPotential() const;
  bool GetAutoPlotKCurve() const;
  
  const char * GetKCurveFile() const;
  const char * GetKFoundFile() const;
  const char * GetPotentialFile() const;

  Potential * GetPotential() const;

  ParametrizedCurve * GetKCurve() const;

  const vector<BasisFunction> & GetBasisFunctions() const;

  ExpectedMatrixType GetExpectedMatrixType() const;

  ///Setters

  void SetAutoPlotPotential(bool value);
  void SetAutoPlotKCurve(bool value);
  
  void SetKCurveFile(const char * value);
  void SetKFoundFile(const char * value);
  void SetPotentialFile(const char * value);
  
  void SetPotential(Potential * value);
  void SetKCurve(ParametrizedCurve * value);
  
  void SetBasisFunctions(vector<BasisFunction> value);




  
  void SetVerbosityLevel(unsigned int value);

  void SetExpectedMatrixType(ExpectedMatrixType val);


private:

  ///This stuff should typically be contained in the configuration file.

  unsigned int verbosityLevel; ///The verbosity level. Default: 0
  bool autoPlotPotential; ///Auto plot potential. Default: false
  bool autoPlotKCurve; ///Auto plot K curve .Default: false
  char kCurveFile[MAX_FILENAME_SIZE]; ///Output file for the basis state curve. Default: KCurve.dat
  char kFoundFile[MAX_FILENAME_SIZE]; ///Output file for found k-values. Default: KFound.dat
  char potentialFile[MAX_FILENAME_SIZE]; ///Output file for potential. Default: Potential.dat

  Potential * potential; ///Potential function. Default: a piecewise potential with 3 nonzero regions.
  ParametrizedCurve * kCurve; ///KCurve. Default: a somewhat standard Berggren curve.
  vector<BasisFunction> basisFunctions; ///Basis functions. Default: two functions, "expi+" and "expi-"

  ExpectedMatrixType matrixType;


  /*
  ///Read out parameters from the command line interpreter (or default ones, if they were not specified.
  int verbosityLevel = atoi((myInterpreter->ReadFlaggedCommandStrict("verbose").front()).c_str());
  unsigned int threads = atoi(myInterpreter->ReadFlaggedCommandStrict("threads").front().c_str());
  double kCutoff = atof(myInterpreter->ReadFlaggedCommandStrict("kCutoff").front().c_str());
  double kMid = atof(myInterpreter->ReadFlaggedCommandStrict("kMid").front().c_str());
  double kDepth = atof(myInterpreter->ReadFlaggedCommandStrict("kDepth").front().c_str());
  unsigned int kValuesFar = atof(myInterpreter->ReadFlaggedCommandStrict("kValuesFar").front().c_str());
  unsigned int kValuesClose = atof(myInterpreter->ReadFlaggedCommandStrict("kValuesClose").front().c_str());
  string dataFile = myInterpreter->ReadFlaggedCommandStrict("dataFile").front().c_str();
  string potentialFile = myInterpreter->ReadFlaggedCommandStrict("potentialFile").front().c_str();


  */

};



#endif
