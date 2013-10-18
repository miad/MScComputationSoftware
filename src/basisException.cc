#include "basisException.hh"

basisException::basisException(string m)
  :msg(m)
{

}

basisException::~basisException() throw()
{

}

const char* basisException::what() const throw()
{
  return msg.c_str();
}
