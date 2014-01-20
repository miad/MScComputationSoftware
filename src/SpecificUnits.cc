#include "SpecificUnits.hh"

SpecificUnits::SpecificUnits(double _hbarTimesLambda, double _massOverLambda2, string _lengthUnitName, string _energyUnitName, double _timeToHertzFactor)
  :hbarTimesLambda(_hbarTimesLambda), massOverLambda2(_massOverLambda2), lengthUnitName(_lengthUnitName), energyUnitName(_energyUnitName), timeToHertzFactor(_timeToHertzFactor)
{
}

SpecificUnits::SpecificUnits()
  :hbarTimesLambda(197.326971812), massOverLambda2(938.0), lengthUnitName("fm"), energyUnitName("MeV"), timeToHertzFactor(1.0)
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


double SpecificUnits::GetTimeToHertzFactor() const
{
  return timeToHertzFactor;
}



ComplexDouble SpecificUnits::EnergyToKValue(const ComplexDouble & energy) const
{
  ///If numerical stability is mean to us, then rotate.
  ComplexDouble toReturn = sqrt(energy * 2. * massOverLambda2) / hbarTimesLambda;
  
  if( (abs(imag(toReturn)) > 1E2*abs(real(toReturn))  && imag(toReturn) < 0))
	{
	  toReturn *= -1.0;
	}

  return toReturn;
}

ComplexDouble SpecificUnits::KValueToEnergy(const ComplexDouble & kValue) const
{
  return pow(kValue * hbarTimesLambda, 2.0) / (2.0 * massOverLambda2);
}

