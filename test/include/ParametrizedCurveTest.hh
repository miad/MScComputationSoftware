#ifndef ParametrizedCurveTest_hh
#define ParametrizedCurveTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "ParametrizedCurve.hh"

#include <iostream>

#ifndef EPSILON
#define EPSILON 1E-9
#endif

#define DBL_EQUAL(d1, d2) (abs((d1)-(d2)) < EPSILON)

using namespace std;

class ParametrizedCurveTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  bool TestCase1() const; /// Test case.
  bool TestCase2() const; /// Test case.
  bool TestCase3() const; /// Test case.
  bool TestCase4() const; /// Test case.
  bool TestCase5() const; /// Test case.
  int TestCase6() const; /// Test case. Returns an error code indicating pos of failure, if fail. Returns 0 if success.
  int TestCase7() const; ///Test case.
  int TestCase8() const; ///Test case.
};
#endif
