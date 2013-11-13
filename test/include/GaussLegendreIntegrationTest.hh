#ifndef GaussLegendreIntegrationTest_hh
#define GaussLegendreIntegrationTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <algorithm>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "LegendreRule.hh"
#include "ParametrizedCurve.hh"

#include <iostream>

#define NPPUSH(v) numberOfPointsOnCurve.push_back(v)


using namespace std;

class GaussLegendreIntegrationTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const;
  int TestCase3() const;
};
#endif
