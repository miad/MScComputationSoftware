#include "InteractionProperties.hh"

InteractionProperties::InteractionProperties()
  :couplingCoefficient(0), precision(1000), nmax(30), cacheFile("")
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

ulong InteractionProperties::GetPrecision() const
{
  return precision;
}

void InteractionProperties::SetPrecision(ulong value)
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

string InteractionProperties::GetCacheFile() const
{
  return cacheFile;
}

void InteractionProperties::SetCacheFile(string value)
{
  cacheFile = value;
}
