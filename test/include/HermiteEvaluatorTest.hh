#ifndef HermiteEvaluatorTest_hh
#define HermiteEvaluatorTest_hh 1

#include <iostream>
#include <cassert>
#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iomanip>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "Matrix.hpp"
#include "HermiteEvaluator.hh"
#include "RLMacros.hpp"
#include "RLException.hh"

#include <algorithm>

using namespace std;

class HermiteEvaluatorTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.

};
#endif
