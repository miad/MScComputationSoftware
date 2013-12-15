#include "SpecificUnits.hh"

SpecificUnits::SpecificUnits(double _hbarTimesLambda, double _massOverLambda2, string _lengthUnitName, string _energyUnitName)
  :hbarTimesLambda(_hbarTimesLambda), massOverLambda2(_massOverLambda2), lengthUnitName(_lengthUnitName), energyUnitName(_energyUnitName)
{
}

SpecificUnits::SpecificUnits()
  :hbarTimesLambda(197.326971812), massOverLambda2(938.0), lengthUnitName("fm"), energyUnitName("MeV")
{
}

SpecificUnits::~SpecificUnits()
{

}

double SpecificUnits::GetHbarTimesLambda() const
{
  return hbarTimesLambda;
}

double SpecificUnits::GetMassOverLambda2() const
{
  return massOverLambda2;
}

string SpecificUnits::GetLengthUnitName() const
{
  return lengthUnitName;
}

string SpecificUnits::GetEnergyUnitName() const
{
  return energyUnitName;
}
