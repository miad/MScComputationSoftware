#ifndef EigenPairTest_hh
#define EigenPairTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "EigenPair.hh"

#include <algorithm>

using namespace std;

class EigenPairTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const;
  int TestCase3() const;
  int TestCase4() const;
  int TestCase5() const;
  int TestCase6() const;
  int TestCase7() const;
  int TestCase8() const;
  int TestCase9() const;

};
#endif
