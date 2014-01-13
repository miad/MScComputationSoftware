#ifndef HarmonicBasisFunctionTest_hh
#define HarmonicBasisFunctionTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "Matrix.hpp"
#include "HarmonicBasisFunction.hh"
#include "RLMacros.hpp"
#include "HermiteRule.hh"
#include "LegendreRule.hh"
#include "ParametrizedPotential.hh"

#include <algorithm>

#define NMAX_CERT_ONE 50
#define NMAX_CERT 50
#define NMAX_CERT_TWO 50

#define KRONECKER_DELTA(n, m) (n==m)
#define GL_EQ(a, b) (abs(a-b)<1E-9)

using namespace std;

class HarmonicBasisFunctionTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const;
  int TestCase3() const;
  int TestCase4() const;
  int TestCase5() const;

};
#endif
