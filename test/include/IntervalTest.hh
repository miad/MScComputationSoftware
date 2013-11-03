#ifndef IntervalTest_hh
#define IntervalTest_hh 1

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
#include "Potential.hh"

#include <algorithm>

#ifndef EPSILON
#define EPSILON 1E-9
#endif

#ifndef PI
#define PI 3.141592653589793238462643
#endif

#define DBL_EQUAL(d1, d2) (abs((d1)-(d2)) < EPSILON)

using namespace std;

class IntervalTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.

};
#endif
