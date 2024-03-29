#ifndef LapackeEigenvalueSolverTest_hh
#define LapackeEigenvalueSolverTest_hh 1

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
#include "LapackeEigenvalueSolver.hh"
#include "EigenPair.hh"

#include <algorithm>

using namespace std;

#define EIGEN_EQUAL(a, b) (DBL_EQUAL(a, b) || DBL_EQUAL(a, -1.0*b))

class LapackeEigenvalueSolverTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const; /// Test case.
  int TestCase3() const;

};
#endif
