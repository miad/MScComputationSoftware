#ifndef SpecificUnits_hh
#define SpecificUnits_hh 1

#include <string>
#include <iostream>

#include "RLMacros.hpp"

using namespace std;

/**
   Class containing the specific units used in the computation process.
 */
class SpecificUnits
{
public:
  SpecificUnits(double _hbarTimesLambda, ///Value
				double _massOverLambda2,  ///Value
				string _lengthUnitName, ///Value 
				string _energyUnitName, ///Value
				double timeToHertzFactor ///Value
				); ///Constructor.
  
  SpecificUnits(); ///Default constructor. 
  ~SpecificUnits();

  double GetHbarTimesLambda() const; /// Getter.
  double GetMassOverLambda2() const; /// Getter. 
  string GetLengthUnitName() const; /// Getter.
  string GetEnergyUnitName() const; /// Getter.

  double GetTimeToHertzFactor() const;
  

  ComplexDouble EnergyToKValue(const ComplexDouble & energy ///The energy.
							   ) const; ///Converts an energy to a wave vector using the supplied internal unit system.

  ComplexDouble KValueToEnergy(const ComplexDouble & kValue ///The value to convert.
							   ) const; ///Converts a wave vector to an energy using the supplied internal unit system.




private:
  double hbarTimesLambda; ///The value of hbar times a constant, lambda.
  double massOverLambda2; ///The value of the mass, divided by lambda (same as hbar) squared.

  string lengthUnitName; ///Name of length unit.
  string energyUnitName; ///Name of energy unit.

  double timeToHertzFactor;
};

#endif
