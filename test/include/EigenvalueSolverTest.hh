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

using namespace std;

#define EIGEN_EQUAL(a, b) (DBL_EQUAL(a, b) || DBL_EQUAL(a, -1.0*b))

class EigenvalueSolverTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const; /// Test case.
  int TestCase3() const;


  static bool CplexCompare(const ComplexDouble & c1,  ///First.
						   const ComplexDouble & c2 ///Second
						   ); ///Compare two complex doubles for sorting. Compares first on real part, then on imaginary part.

};
#endif
