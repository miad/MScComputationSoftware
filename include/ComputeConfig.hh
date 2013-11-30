#ifndef ComputeConfig_hh
#define ComputeConfig_hh 1

#include <list>
#include <vector>


#include "libconfig.hh"
#include "RLException.hh"
#include "Globals.hpp"
#include "Potential.hh"
#include "ParametrizedCurve.hh"
#include "BasisFunction.hh"
#include "RLMacros.hpp"
#include "PiecewiseConstantPotential.hh"


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

  unsigned int GetNumberOfThreads() const;

  bool GetAutoPlotPotential() const;
  bool GetAutoPlotKCurve() const;
  
  const char * GetKCurveFile() const;
  const char * GetKFoundFile() const;
  const char * GetPotentialFile() const;
  const char * GetPotentialPrecisionFile() const;

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
  void SetPotentialPrecisionFile(const char * value);
  
  void SetPotential(Potential * value);
  void SetKCurve(ParametrizedCurve * value);
  
  void SetBasisFunctions(vector<BasisFunction> value);




  
  void SetVerbosityLevel(unsigned int value);
  void SetNumberOfThreads(unsigned int value);

  void SetExpectedMatrixType(ExpectedMatrixType val);


private:

  ///This stuff should typically be contained in the configuration file.
  
  unsigned int numberOfThreads; ///number of threads. Default: 1
  unsigned int verbosityLevel; ///The verbosity level. Default: 0
  bool autoPlotPotential; ///Auto plot potential. Default: false
  bool autoPlotKCurve; ///Auto plot K curve .Default: false
  char kCurveFile[MAX_FILENAME_SIZE]; ///Output file for the basis state curve. Default: KCurve.dat
  char kFoundFile[MAX_FILENAME_SIZE]; ///Output file for found k-values. Default: KFound.dat
  char potentialFile[MAX_FILENAME_SIZE]; ///Output file for potential. Default: Potential.dat
  char potentialPrecisionFile[MAX_FILENAME_SIZE]; ///Output file for potential precision points. Default: PotentialPrec.dat

  Potential * potential; ///Potential function. Default: a piecewise potential with 3 nonzero regions.
  ParametrizedCurve * kCurve; ///KCurve. Default: a somewhat standard Berggren curve.
  vector<BasisFunction> basisFunctions; ///Basis functions. Default: two functions, "expi+" and "expi-"

  ExpectedMatrixType matrixType;
};

#endif
