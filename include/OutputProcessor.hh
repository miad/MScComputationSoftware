#ifndef OutputProcessor_hh
#define OutputProcessor_hh 1

#include <algorithm>
#include <cmath>
#include <iostream>

#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "OutputFilenames.hh"
#include "LapackeEigenvalueSolver.hh"
#include "RLException.hh"
#include "ComputeConfig.hh"
#include "EigenInformation.hh"
#include "Matrix.hpp"
#include "MultiTasker.hpp"
#include "WavefunctionTwoParticleWorkerData.hh"
#include "RLMacros.hpp"
#include "CompositeBasisFunction.hh"



/** This class will process the results from computations, and save it in the correct format in the correct output files. In order to do so it has to do some filtering etc. of points.

 */
class OutputProcessor : public VerbosePrinterEventEnabled
{
public:
  OutputProcessor(ComputeConfig * _config ///Ownership is NOT passed. The object may NOT be destroyed before this class to guarantee integrity.
				  ); ///Constructor.

  ~OutputProcessor(); ///Destructor.

  void SetEigenInformation(EigenInformation * data ///Eigeninfo to set.
						   ); ///Set eigeninfo. Must be done before calling WriteOutput. NOTE: Ownership is not passed, and data may not be deleted before this object is destroyed.
  
  void WritePostOutput() const; ///Call this to perform all POST output processing on the set value of data. Throws exception if data is not set.

  void SaveMatrix(CMatrix * toSave ///The matrix.
				  ) const; ///Write the matrix to predetermined file.

  void SetCompositeBasisFunctions(vector<CompositeBasisFunction * > * value
								  );



private:

  static bool ComplexCompare(const ComplexDouble & d1, const ComplexDouble & d2);



  // Performs certain portions of the output writing, etc.
  void WritePotentialToFile() const; ///Writes output.
  void WritePotentialPrecisionToFile() const; ///Writes output. 
  void WriteKCurveToFile() const; ///Writes output.
  void WriteKFoundToFile() const; ///Writes output.
  void WriteEnergiesToFile() const; ///Writes output.
  void WriteInterestingKPointsVerbosely() const; ///Writes output.
  void WriteInterestingKPointsToFile() const; ///Writes output.
  void WriteInterestingOneParticleWavefunctionsToFile() const; ///Writes output.


  void WriteInterestingTwoParticleWavefunctionsToFile() const;
  void WriteInterestingRelativeTwoParticleWavefunctionsToFile() const;
  void WriteProductTwoParticleWavefunctionToFile() const;


  // Auxiliary functions.

  vector<uint> FindInterestingKPointIndex() const; ///Returns the interesting K points by index (in the eigenData object)

  vector<ComplexDouble> FindInterestingKPoints() const; ///Returns the interesting K points by value.


  static FILE * AssuredFopen(const string filename ///File to open.
							 ); ///Opens a file FOR WRITING and returns the handle, throws exception if the open failed. Ownership passed on.


  vector<double> GetOneParticleBasisRatio(uint eigenIndex, ///Index of wavevector to compute ratio for.
										  double & totalSum /// Output value. May be used.
										  ) const; ///Returns a vector with entries corresponding to the relative dominance (all sum up to 1) of each basis vector


  static double SquareSum(const vector<ComplexDouble> & toSum 
						  );



  
private:
  ComputeConfig * config;
  EigenInformation * eigenData;
  vector<CompositeBasisFunction * > * myCompositeBasisFunctions;

};

#endif
