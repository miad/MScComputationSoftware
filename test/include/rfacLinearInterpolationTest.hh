#ifndef rfacLinearInterpolationTest_hh
#define rfacLinearInterpolationTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>

#include "GenericUnitTest.hh"
#include "rfacLinearInterpolation.hh"
#include "rfacException.hh"

#ifndef EPSILON
#define EPSILON 1E-9
#endif

using namespace std;

class rfacLinearInterpolationTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  bool rfacLinearInterpolation_ConstructsCorrectly_AssertTrue() const; /// Test case.
  bool rfacLinearInterpolation_ThrowsOnTooSmallInput_AssertTrue() const; /// Test case.
  bool rfacLinearInterpolation_PerformsCorrectFitI_AssertTrue(double i) const;

  rfacLinearInterpolation initializeDefault1() const;
};
#endif
