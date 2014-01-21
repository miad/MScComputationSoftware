#ifndef InteractionProperties_hh
#define InteractionProperties_hh 1

#include "RLMacros.hpp"

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

  long GetPrecision() const;

  void SetPrecision(long value
					);

private:
  double couplingCoefficient;
  long precision;
  uint nmax;
};

#endif
