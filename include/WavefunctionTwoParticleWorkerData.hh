#ifndef WavefunctionTwoParticleWorkerData_hh
#define WavefunctionTwoParticleWorkerData_hh 1

#include <vector>
#include "RLMacros.hpp"
#include "ComputeConfig.hh"

/*
  Utility class for worker data when computing 2d wavefunctions (heavy job) in parallell.
  
 */
class WavefunctionTwoParticleWorkerData
{
public:
  WavefunctionTwoParticleWorkerData(vector<ComplexDouble> * _eigVect,
									vector<vector<ComplexDouble> > * _wavefunctionValues,
									ComputeConfig * _config,
									uint _eigenIndex
									);
  ~WavefunctionTwoParticleWorkerData();

  vector<ComplexDouble> * eigVect;
  vector<vector<ComplexDouble> > * wavefunctionValues;
  ComputeConfig * config;
  uint eigenIndex;
private:
  ///Prevent copying, since it with current implementation causes glibc.
  WavefunctionTwoParticleWorkerData(const WavefunctionTwoParticleWorkerData &);
  WavefunctionTwoParticleWorkerData & operator=(const WavefunctionTwoParticleWorkerData &);
};


#endif

