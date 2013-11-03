#ifndef MatrixTest_hh
#define MatrixTest_hh 1

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

#include <iostream>

#ifndef EPSILON
#define EPSILON 1E-9
#endif

#define DBL_EQUAL(d1, d2) (abs((d1)-(d2)) < EPSILON)

using namespace std;

class MatrixTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const; ///Test case.
  int TestCase3() const;
};
#endif
