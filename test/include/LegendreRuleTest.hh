#ifndef LegendreRuleTest_hh
#define LegendreRuleTest_hh 1

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

#include <iostream>
using namespace std;

class LegendreRuleTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  bool TestCase1() const; /// Test case.
  bool TestCase2() const; /// Test case.
  bool TestCase3() const; /// Test case.
  bool TestCase4() const; /// Test case.
};
#endif
