#ifndef rfacDecayableObjectTest_hh
#define rfacDecayableObjectTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>

#include "GenericUnitTest.hh"
#include "rfacDecayableObject.hh"
#include "rfacException.hh"
#include "QuantumNumbers.hh"

#include <iostream>

#ifndef EPSILON
#define EPSILON 1E-9
#endif

using namespace std;

class rfacDecayableObjectTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  bool rfacDecayableObject_ConstructsCorrectly_AssertTrue() const; /// Test case.
  bool rfacDecayableObject_DecaysWithoutException_AssertTrue() const; /// Test case.

};
#endif
