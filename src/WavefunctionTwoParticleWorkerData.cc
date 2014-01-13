#include "WavefunctionTwoParticleWorkerData.hh"


WavefunctionTwoParticleWorkerData::WavefunctionTwoParticleWorkerData(
																	 vector<ComplexDouble> * _eigVect,
																	 vector<vector<ComplexDouble> > * _wavefunctionValues,
																	 ComputeConfig * _config,
																	 uint _eigenIndex
								  )
  :eigVect(_eigVect), wavefunctionValues(_wavefunctionValues), config(_config), eigenIndex(_eigenIndex)
{

}

WavefunctionTwoParticleWorkerData::~WavefunctionTwoParticleWorkerData()
{
  /// Ownership was passed, so delete it.
  delete eigVect;
  eigVect = NULL;
}
