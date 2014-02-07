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
#include "HarmonicBasisFunction.hh"
#include "InteractionProperties.hh"
#include "EigenvalueSolver.hh"


using namespace std;
using namespace libconfig;

#define CONFIG_FILE_VERSION 0.40


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

  uint GetVerbosityLevel() const;

  uint GetNumberOfThreads() const;

  bool GetAutoPlotPotential() const;
  bool GetAutoPlotKCurve() const;
  bool GetAutoPlotWavefunctions() const;

  const OutputFilenames * GetOutputFilenames() const;
  const SpecificUnits * GetSpecificUnits() const;
  

  Potential * GetPotential(uint index = 0 ///Index of particle.
						   ) const;

  ParametrizedCurve * GetKCurve() const;

  const vector<BasisFunction> & GetBasisFunctions() const;

  const vector<ComplexDouble> & GetExtraInterestingPoints() const;


  HarmonicBasisFunction * GetHarmonicBasisFunction() const;

  void SetHarmonicBasisFunction(HarmonicBasisFunction * value);

  ExpectedMatrixType GetExpectedMatrixType() const;


  double GetMinWavefunctionX() const;
  double GetMaxWavefunctionX() const;
  double GetWavefunctionStepsizeX() const;

  uint GetNumberOfParticles() const;

  const InteractionProperties * GetInteractionProperties() const;

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
  
  void SetPotential(Potential * value,
					uint index ///Index of particle.
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


  
  void SetVerbosityLevel(uint value
						 );
  void SetNumberOfThreads(uint value
						  );

  void SetExpectedMatrixType(ExpectedMatrixType val
							 );

  void SetNumberOfParticles(uint value
							);
  
  void SetInteractionProperties(InteractionProperties value
								);

  void SetCouplingCoefficient(double value
							  );

  bool GetHarmonicOverride() const;
  uint GetHarmonicNmax() const;
  void SetHarmonicOverride(bool value
						   );


  void SetHarmonicNmax(uint value
					   );


  EigenvalueSolver * GetSolver();

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

  void ReadMultiParticleData(Setting & computation
							 );
  
  void ReadHarmonicOscillator(Setting & computation
							  );


  void ReadPiecewiseConstantPotential(Setting & poten
									  );

  void ReadParametrizedPotential(Setting & poten
								 );

  void ReadInteractionProperties(Setting & computation
								 );

  void ReadSolverInfo(Setting & computation
					  );


private:

  vector<ComplexDouble> extraInterestingPoints; ///An extra filtering option for artificially inducing extra InterestingPoints to the output, main use is quality control purposes.

  ///This stuff should typically be contained in the configuration file.

  uint numberOfThreads; ///number of threads. Default: 1
  uint verbosityLevel; ///The verbosity level. Default: 0

  bool autoPlotPotential; ///Auto plot potential. Default: false
  bool autoPlotKCurve; ///Auto plot K curve .Default: false
  bool autoPlotWavefunctions; ///Auto plot wavefunctions. Default: false

  double minWavefunctionX; ///Min wavefunction x.
  double maxWavefunctionX; ///Max wavefunction x.
  double wavefunctionStepsizeX; ///Wavefunction stepsize x.

  uint numberOfParticles; ///Number of particles in many-particle case.


  bool harmonicOverride; ///Harmonic override (use harmonic basis instead)
  uint harmonicNmax; ///Use only in case of override.
 

  EigenvalueSolver mySolver;

  vector<Potential*> potentials; ///One for each particle...

  HarmonicBasisFunction * myHarmonicBasisFunction;

  ParametrizedCurve * kCurve; ///KCurve. Default: a somewhat standard Berggren curve.
  vector<BasisFunction> basisFunctions; ///Basis functions. 

  OutputFilenames outputFilenames; ///The output filenames.
  SpecificUnits specificUnits; /// Units to use in computation.

  ExpectedMatrixType matrixType; ///The expected matrix type.

  InteractionProperties myInteractionProperties;
};

#endif
