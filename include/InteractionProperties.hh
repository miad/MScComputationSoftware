#ifndef InteractionProperties_hh
#define InteractionProperties_hh 1

#include "RLMacros.hpp"
#include <string>
using namespace std;

/** Container class for some interaction term properties in 2-body Hamiltonian.
 */
class InteractionProperties
{
public:
  InteractionProperties();

  ~InteractionProperties();

  double GetCouplingCoefficient() const;

  void SetCouplingCoefficient(double value
							  );

  uint GetNMax() const;

  void SetNMax(uint value
			   );

  double GetLowerIntegrationLimit() const;

  void SetLowerIntegrationLimit(double value
								);

  double GetUpperIntegrationLimit() const;

  void SetUpperIntegrationLimit(double value
								);

  ulong GetPrecision() const;

  void SetPrecision(ulong value
					);

  string GetCacheFile() const;

  void SetCacheFile(string value
					);

private:
  double couplingCoefficient;
  ulong precision;
  uint nmax;
  string cacheFile;
};

#endif
