#ifndef PiecewiseConstantPotentialTest_hh
#define PiecewiseConstantPotentialTest_hh 1

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
#include "PiecewiseConstantPotential.hh"

#include <algorithm>

using namespace std;

class PiecewiseConstantPotentialTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const; /// Test case.
  int TestCase3() const; /// Test case.

};
#endif
