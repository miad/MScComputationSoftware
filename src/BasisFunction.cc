#include "BasisFunction.hh"

BasisFunction::BasisFunction(const char * _name)
  :type(-1)
{
  factor = ComplexDouble(1, 0);
  unsigned short len = strlen(_name);
  if(len < 3 || len > 5)
	throw RLException("Invalid name length: %d, name: '%s'", len, _name);

  strcpy(name, _name);


  for(int i = 0; i<len; ++i)
	name[i] = tolower(name[i]);

  char fun[4];
  memset(fun, '\0', sizeof(char)*4);

  strncpy(fun, name, 3);
  if(strcmp(fun, "exp") == 0)
	type = 0;
  if(strcmp(fun, "sin") == 0)
	type = 1;
  if(strcmp(fun, "cos") == 0)
	type = 2;
  if(type == -1)
	throw RLException("Unknown function type: %s", fun);

  for(int j = 3; j<len; ++j)
	{
	  if(name[j] == '+' || name[j] == '\n' || name[j] == ' ')
		{
		  ///Acceptable, but no action taken.
		}
	  else if(name[j] == '-')
		{
		  factor *= -1.0;
		}
	  else if(name[j] == 'i')
		{
		  factor *= ComplexDouble(0, 1);
		}
	  else
		{
		  throw RLException("Unknown stuff in the name string sent to BasisFunction constructor.");
		}
	}
}

ComplexDouble BasisFunction::Eval(ComplexDouble x)
{
  x *= factor;
  switch(type)
	{
	case 0:
	  return exp(x);
	case 1:
	  return sin(x);
	case 2:
	  return cos(x);
	default:
	  throw RLException("Invalid type. This is an internal inconsistency in the code for BasisFunction.");
	}
}

const char * BasisFunction::GetName() const
{
  return name;
}

double BasisFunction::GetPreFactor() const
{
  switch(type)
	{
	case 0:
	  return 1;
	case 1:
	  return 2;
	case 2:
	  return 2;
	default:
	  throw RLException("Invalid type. This is an internal inconsistency in the code for BasisFunction.");
	}
}
