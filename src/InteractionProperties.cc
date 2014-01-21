#include "InteractionProperties.hh"

InteractionProperties::InteractionProperties()
  :couplingCoefficient(0), precision(1000), nmax(30)
{

}

InteractionProperties::~InteractionProperties()
{

}


double InteractionProperties::GetCouplingCoefficient() const
{
  return couplingCoefficient;
}

void InteractionProperties::SetCouplingCoefficient(double value)
{
  couplingCoefficient = value;
}

long InteractionProperties::GetPrecision() const
{
  return precision;
}

void InteractionProperties::SetPrecision(long value)
{
  precision = value;
}

uint InteractionProperties::GetNMax() const
{
  return nmax;
}

void InteractionProperties::SetNMax(uint value)
{
  nmax = value;
}
