#include "Potential.hh"

Potential::Potential()
  :minX(0.), maxX(0.)
{

}

Potential::Potential(string filename)
{
  FILE * inputFile = fopen(filename.c_str(), "r");
  char * line = NULL; ///stores line data when read from file.
  size_t len = 0; ///is ignored. 
  ssize_t numberOfReadCharacters; ///length of read text

  if(inputFile == NULL)
	throw RLException("Could not open potential input file '%s'.", filename.c_str());

  while ((numberOfReadCharacters = getline(&line, &len, inputFile)) != -1) 
	{
	  ///check if we should ignore the line.
	  int ptr = 0;
	  while(ptr < numberOfReadCharacters && line[ptr] == ' ')
		++ptr;
	  if(ptr == numberOfReadCharacters || line[ptr] == '#' || line[ptr] == '\n')
		continue;

	  double x1, x2, y;
	  if(sscanf(line, "%lf %lf %lf",&x1, &x2, &y) == EOF)
		{
		  throw RLException("Suspicious line in potential file '%s': '%s'", filename.c_str(), line);
		}
	  AddValue(x1, x2, y);
	}

  if (line)
	free(line);
  fclose(inputFile);
  if(PotentialPoints.empty())
	throw RLException("The specified potential file '%s' does not contain any potential data.", filename.c_str());
}

Potential::~Potential()
{

}

ComplexDouble Potential::FastExpIntegrate(ComplexDouble exponentVal) const
{
  ComplexDouble value = 0;
  if(abs(exponentVal) < EPS) ///If the exponent is zero, the integrand is computed differently.
	{
	  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
		{
		  value += it->y*(it->x2 - it->x1);
		}
	  return value;
	}

  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{
	  value+=it->y*(exp(ComplexDouble(0, 1)*exponentVal*it->x2) - 
					exp(ComplexDouble(0, 1)*exponentVal*it->x1) );
	}
  value /= ComplexDouble(0, 1)*exponentVal;

  return value;
}

ComplexDouble Potential::FastCosIntegrate(ComplexDouble k1, ComplexDouble k2) const
{  
  ComplexDouble value = 0;
  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{
	  if( DBL_EQUAL(k1, k2) && DBL_EQUAL(k1, 0.)) ///If the exponent is zero, the integrand is computed differently.
		{
		  value += it->y*(it->x2 - it->x1);
		}
	  else if( DBL_EQUAL(k1, k2) || DBL_EQUAL(k1, -k2) )
		{
		  value += it->y * ( (it->x2 - it->x1)*0.5 +
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else
		{
		  value += it->y * 0.5*(
								(sin((k1-k2)*it->x2)-sin((k1-k2)*it->x1))/(k1-k2) +
								(sin((k1+k2)*it->x2)-sin((k1+k2)*it->x1))/(k1+k2)
								);
		}
	}
  return value;
}

ComplexDouble Potential::FastSinIntegrate(ComplexDouble k1, ComplexDouble k2) const
{

  ComplexDouble value = 0;
  for(list<Interval>::const_iterator it = PotentialPoints.begin(); it!=PotentialPoints.end(); ++it)
	{	  
	  if( DBL_EQUAL(k1, k2) && DBL_EQUAL(k1, 0.)) ///If the exponent is zero, the integrand is computed differently.
		{
		  value = 0.;
		  return value;
		}
	  else if( DBL_EQUAL(k1, k2))
		{

		  value += it->y * ( (it->x2 - it->x1)*0.5 -
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else if( DBL_EQUAL(k1, -k2))
		{

		  value += -1.*it->y * ( (it->x2 - it->x1)*0.5 -
								   (sin(2.*k1*it->x2)-sin(2.*k1*it->x1)) /(4.*k1)
								  );
		}
	  else
		{

		  value += it->y * 0.5*(
								(sin((k1-k2)*it->x2)-sin((k1-k2)*it->x1))/(k1-k2) -
								(sin((k1+k2)*it->x2)-sin((k1+k2)*it->x1))/(k1+k2)
								);
		}
	}
  return value;
}


void Potential::AddValue(Interval toAdd)
{
  foric(list<Interval>, PotentialPoints, it)
	{
	  if(it->Overlaps(toAdd))
		{
		  throw RLException("Attempted to add an overlapping interval.");
		}
	}
  minX = PotentialPoints.empty()?toAdd.x1:MIN(minX, toAdd.x1);
  maxX = PotentialPoints.empty()?toAdd.x2:MAX(maxX, toAdd.x2);

  PotentialPoints.push_back(toAdd);
}

void Potential::AddValue(double x1, double x2, double y)
{
  AddValue(Interval(x1, x2, y));
}

double Potential::Evaluate(double x) const
{
  foric(list<Interval>, PotentialPoints, it)
	{
	  if( it->x1 <= x && it->x2 > x )
		{
		  return it->y;
		}
	}
  return 0;
}

list<Interval> Potential::GetPotentialPoints() const
{
  return PotentialPoints;
}

double Potential::GetMinX() const
{
  return minX;
}

double Potential::GetMaxX() const
{
  return maxX;
}

unsigned int Potential::GetNumberOfValues() const
{
  return PotentialPoints.size();
}


void Potential::Clear() 
{
  PotentialPoints.clear();
  minX = 1E99;
  maxX = -1E99;
}


void Potential::BasisIntegrate(BasisFunction & b1, BasisFunction & b2, ComplexDouble & k1, ComplexDouble & k2)
{
  
}
