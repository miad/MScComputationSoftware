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
#include "ParametrizedPotential.hh"
#include "OutputFilenames.hh"
#include "SpecificUnits.hh"


using namespace std;
using namespace libconfig;

#define CONFIG_FILE_VERSION 0.35


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
  bool GetAutoPlotWavefunctions() const;

  const OutputFilenames * GetOutputFilenames() const;
  const SpecificUnits * GetSpecificUnits() const;
  

  Potential * GetPotential() const;

  ParametrizedCurve * GetKCurve() const;

  const vector<BasisFunction> & GetBasisFunctions() const;

  const vector<ComplexDouble> & GetExtraInterestingPoints() const;

  ExpectedMatrixType GetExpectedMatrixType() const;


  double GetMinWavefunctionX() const;
  double GetMaxWavefunctionX() const;
  double GetWavefunctionStepsizeX() const;


  ///Setters


  void SetExtraInterestingPoints(vector<ComplexDouble> value
								 );

  void SetAutoPlotPotential(bool value
							);
  void SetAutoPlotKCurve(bool value
						 );
  void SetAutoPlotWavefunctions(bool value
								);


  void SetOutputFilenames(const OutputFilenames & value
						  );

  void SetSpecificUnits(const SpecificUnits & value
						);
  
  void SetPotential(Potential * value
					);
  void SetKCurve(ParametrizedCurve * value
				 );
  
  void SetBasisFunctions(vector<BasisFunction> value
						 );

  void SetMinWavefunctionX(double value
						   );
  void SetMaxWavefunctionX(double value
						   );
  void SetWavefunctionStepsizeX(double value
								);


  
  void SetVerbosityLevel(unsigned int value
						 );
  void SetNumberOfThreads(unsigned int value
						  );

  void SetExpectedMatrixType(ExpectedMatrixType val
							 );


private:
  void ReadOutputFiles(Setting & output
					   );
  void ReadExpectedMatrixType(Setting & computation
							  );
  void ReadSpecificUnits(Setting & computation
						 );
  void ReadPotential(Setting & computation
					 );
  void ReadKCurve(Setting & computation
				  );
  void ReadBasisFunctions(Setting & computation
						  );

  void ReadAutoLaunch(Setting & program
					  );

  void ReadWavefunctionProperties(Setting & outputSpecifics
								  );

  void ReadExtraInteresting(Setting & outputSpecifics
							);


private:

  vector<ComplexDouble> extraInterestingPoints; ///An extra filtering option for artificially inducing extra InterestingPoints to the output, main use is quality control purposes.

  ///This stuff should typically be contained in the configuration file.

  unsigned int numberOfThreads; ///number of threads. Default: 1
  unsigned int verbosityLevel; ///The verbosity level. Default: 0

  bool autoPlotPotential; ///Auto plot potential. Default: false
  bool autoPlotKCurve; ///Auto plot K curve .Default: false
  bool autoPlotWavefunctions; ///Auto plot wavefunctions. Default: false

  double minWavefunctionX;
  double maxWavefunctionX;
  double wavefunctionStepsizeX;


  Potential * potential; ///Potential function. Default: a piecewise potential with 3 nonzero regions.
  ParametrizedCurve * kCurve; ///KCurve. Default: a somewhat standard Berggren curve.
  vector<BasisFunction> basisFunctions; ///Basis functions. Default: two functions, "expi+" and "expi-"

  OutputFilenames outputFilenames; ///The output filenames.
  SpecificUnits specificUnits; /// Units to use in computation.

  ExpectedMatrixType matrixType; ///The expected matrix type.
};

#endif
