#ifndef EigenvalueSolverTest_hh
#define EigenvalueSolverTest_hh 1

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
#include "EigenvalueSolver.hh"

#include <algorithm>

#ifndef EPSILON
#define EPSILON 1E-9
#endif

#define DBL_EQUAL(d1, d2) (abs((d1)-(d2)) < EPSILON)

using namespace std;

class EigenvalueSolverTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const; /// Test case.


  static bool CplexCompare(const ComplexDouble & c1,  ///First.
						   const ComplexDouble & c2 ///Second
						   ); ///Compare two complex doubles for sorting. Compares first on real part, then on imaginary part.

};
#endif
