#ifndef OutputProcessor_hh
#define OutputProcessor_hh 1

#include <algorithm>
#include <cmath>
#include <iostream>

#include "VerbosePrinter.hh"
#include "VerbosePrinterEventEnabled.hh"
#include "OutputFilenames.hh"
#include "EigenvalueSolver.hh"
#include "RLException.hh"
#include "ComputeConfig.hh"

/** This class will process the results from computations, and save it in the correct format in the correct output files. In order to do so it has to do some filtering etc. of points.

 */
class OutputProcessor : public VerbosePrinterEventEnabled
{
public:
  OutputProcessor(ComputeConfig * _config ///Ownership is NOT passed. The object may NOT be destroyed before this class to guarantee integrity.
				  ); ///Constructor.

  ~OutputProcessor(); ///Destructor.

  void SetEigenInformation(EigenInformation * data);
  
  void WriteOutput() const; ///Call this to perform ALL output processing on the set value of data. Throws exception if data is not set.



private:

  // Performs certain portions of the output writing, etc.
  void WritePotentialToFile() const; ///Writes output.
  void WritePotentialPrecisionToFile() const; ///Writes output. 
  void WriteKCurveToFile() const; ///Writes output.
  void WriteKFoundToFile() const; ///Writes output.
  void WriteInterestingKPointsVerbosely() const; ///Writes output.
  void WriteInterestingKPointsToFile() const; ///Writes output.
  vector<ComplexDouble> FindInterestingKPoints() const; ///Writes output.
  void WriteInterestingWavefunctionsToFile() const; ///Writes output.


  // Auxiliary functions.

  static FILE * AssuredFopen(const string filename ///File to open.
							 ); ///Opens a file FOR WRITING and returns the handle, throws exception if the open failed. Ownership passed on.

  ComplexDouble EnergyToKValue(const ComplexDouble & energy ///The energy.
							   ) const; ///Converts an energy to a wave vector using the supplied internal unit system.

  ComplexDouble KValueToEnergy(const ComplexDouble & kValue ///The value to convert.
							   ) const; ///Converts a wave vector to an energy using the supplied internal unit system.


  vector<double> GetBasisRatio(unsigned int eigenIndex, ///Index of wavevector to compute ratio for.
							   double & totalSum /// Output value. May be used.
							   ) const; ///Returns a vector with entries corresponding to the relative dominance (all sum up to 1) of each basis vector




  
private:
  ComputeConfig * config;
  EigenInformation * eigenData;


  




};

#endif
