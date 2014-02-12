#ifndef HarmonicBasisFunction_hh
#define HarmonicBasisFunction_hh 1

#include "Globals.hpp"
#include "RLException.hh"
#include "fparser.hh"
#include <stdio.h>
#include <ctype.h>
#include <cstring>
#include <map>
#include <iostream>
#include <vector>
#include <utility>
#include "HermiteEvaluator.hh"
#include "HermiteRule.hh"
#include "Potential.hh"
#include "SpecificUnits.hh"
#include "LegendreRule.hh"
#include "RLMacros.hpp"

#define MAX_FACTORIALS 300

#define TOO_SMALL(x) (abs(x) < 1E-300)
#define TOO_LARGE(x) (abs(x) > 1E300)
#ifdef MASS
#error MASS macro was already defined.
#endif
#define MASS (units->GetMassOverLambda2())
#ifdef HBAR
#error HBAR macro already defined.
#endif
#define HBAR (units->GetHbarTimesLambda())

using namespace std;


class HarmonicBasisFunction
{
public:
  HarmonicBasisFunction(vector<double> _xmin,
						vector<double> _omega,
						vector<Potential *> * _potentials,
						SpecificUnits * _units,
						uint precision = 300
						);
						
  ~HarmonicBasisFunction();


  long double DiffIntegrate(uint n1, 
							uint n2,
							uint pIndex
							) const;

  long double Integrate(uint n1, 
						uint n2,
						uint pIndex
						) const;

  long double NormalizationConstant(uint n, uint pIndex) const;

  long double Eval(uint n, double x, uint pIndex) const;

  long double EvalNonExponentPart(uint n, double x, uint pIndex) const;

  long double KineticTerm(uint n1, 
						  uint n2,
						  uint pIndex
						  ) const;
  
  double GetEigenEnergy(uint n,
						uint pIndex
						) const;

  double GetXmin(uint pIndex) const;
  double GetOmega(uint pIndex) const;

private:

  void InitFactorials();

  vector<long double> factorials;
  vector<vector<pair<double, double> > > GHpoints;
  vector<double> xmin;
  vector<double> omega;
  vector<Potential *> * potentials;
  SpecificUnits * units;
};

#endif
