#ifndef GaussHermiteIntegrationTest_hh
#define GaussHermiteIntegrationTest_hh 1

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
#include "HermiteRule.hh"
#include "HermiteEvaluator.hh"

#include <iostream>

using namespace std;

class GaussHermiteIntegrationTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const;

};
#endif
