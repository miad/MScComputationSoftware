#include "rfacDecayableObjectTest.hh"

bool rfacDecayableObjectTest::rfacDecayableObject_ConstructsCorrectly_AssertTrue() const
{
  rfacDecayableObject * myDecayableObject = new rfacDecayableObject();
  myDecayableObject->SetQuantumNumbers(QuantumNumbers(5,3,5,5,1));
  myDecayableObject->SetStateFile("testsample/sample1.root");
  myDecayableObject->InitializeDecay();
  return true;
}

bool rfacDecayableObjectTest::rfacDecayableObject_DecaysWithoutException_AssertTrue() const
{
  rfacDecayableObject * myDecayableObject = new rfacDecayableObject();
  myDecayableObject->SetQuantumNumbers(QuantumNumbers(7,1,1,1,1));
  myDecayableObject->SetStateFile("testsample/sample1.root");
  myDecayableObject->InitializeDecay();
  myDecayableObject->DoStepTime(1E-5); //since this is stocastic, don't check it.
  return true;
}


int rfacDecayableObjectTest::runUnitTests() const
{
  cout << "Running unit tests on rfacDecayableObject...";
  cout << flush;
  if(!rfacDecayableObject_ConstructsCorrectly_AssertTrue())
    return 1;
  if(!rfacDecayableObject_DecaysWithoutException_AssertTrue())
    return 2;
  cout << "done" << endl;
  return 0;
}
