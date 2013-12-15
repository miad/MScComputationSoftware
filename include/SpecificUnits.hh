#ifndef SpecificUnits_hh
#define SpecificUnits_hh 1

#include <string>

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
				string _energyUnitName ///Value
				); ///Constructor.
  
  SpecificUnits(); ///Default constructor. 

  double GetHbarTimesLambda() const; /// Getter.
  double GetMassOverLambda2() const; /// Getter. 
  string GetLengthUnitName() const; /// Getter.
  string GetEnergyUnitName() const; /// Getter.


private:
  double hbarTimesLambda; ///The value of hbar times a constant, lambda.
  double massOverLambda2; ///The value of the mass, divided by lambda (same as hbar) squared.

  string lengthUnitName; ///Name of length unit.
  string energyUnitName; ///Name of energy unit.
};

#endif
