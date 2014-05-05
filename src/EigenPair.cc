#include "EigenPair.hh"

EigenPair::EigenPair()
{
  Eigenvalue = 0.0;
}

EigenPair::EigenPair(size_t size)
{
  Eigenvalue = 0.0;
  Eigenvector.resize(size, 0.0);
}

int EigenPair::ComplexFirst(const ComplexDouble & c1, const ComplexDouble & c2)
{
  //we need to somehow define an ordering of the eigenpairs. Arbitrarily:
  if(c1.real() < c2.real())
	return 1;
  else if(c1.real() > c2.real() )
	return -1;

  if(c1.imag() < c2.imag())
	return 1;
  else if(c1.imag() > c2.imag() )
	return -1;
  return 0;
}

bool EigenPair::operator<(const EigenPair & other) const
{
  int cfr = ComplexFirst(Eigenvalue, other.Eigenvalue);
  if(cfr)
	return cfr > 0;

  if(Eigenvector.size() < other.Eigenvector.size())
	return true;
  else if(Eigenvector.size() > other.Eigenvector.size())
	return false;
  
  for(size_t i = 0; i<Eigenvector.size(); ++i)
	{
	  int cfr2 = ComplexFirst(Eigenvector[i], other.Eigenvector[i]);
	  if(cfr2)
		return cfr2 > 0;
	}

  return false;
}
