#include "BasisFunctionTest.hh"

int BasisFunctionTest::TestCase1() const
{
  BasisFunction Fsin("sin(k*x)");
  BasisFunction Fsinp("sin(1*k*x)");
  BasisFunction Fsinip("sin(1i*k*x)");
  BasisFunction Fcos("cos(k*x)");
  BasisFunction Fcosp("cos(k*x)");
  BasisFunction Fcosip("cos(1i*k*x)");
  BasisFunction Fexp("exp(x)");
  BasisFunction Fexpp("exp(x)");
  BasisFunction Fexpm("exp(-x)");
  BasisFunction Fexpip("exp(1i*x)");
  BasisFunction Fexpim("exp(-1i*x)");
  
  if(strcmp(Fsin.GetName(), "sin(k*x)") != 0)
	return 2;
  if(strcmp(Fsinp.GetName(), "sin(1*k*x)") != 0)
	return 3;
  if(strcmp(Fsinip.GetName(), "sin(1i*k*x)") != 0)
	return 4;
  if(strcmp(Fcos.GetName(), "cos(k*x)") != 0)
	return 5;
  if(strcmp(Fcosp.GetName(), "cos(k*x)") != 0)
	return 6;
  if(strcmp(Fcosip.GetName(), "cos(1i*k*x)") != 0)
	return 7;
  if(strcmp(Fexp.GetName(), "exp(x)") != 0)
	return 8;
  if(strcmp(Fexpp.GetName(), "exp(x)") != 0)
	return 9;
  if(strcmp(Fexpm.GetName(), "exp(-x)") != 0)
	return 10;
  if(strcmp(Fexpip.GetName(), "exp(1i*x)") != 0)
	return 11;
  if(strcmp(Fexpim.GetName(), "exp(-1i*x)") != 0)
	return 12;

  if(!DBL_EQUAL(Fsin.Eval(1.42), sin(1.42)))
	return 13;
  if(!DBL_EQUAL(Fsinp.Eval(1.42), sin(1.42)))
	return 14;
  if(!DBL_EQUAL(Fsinip.Eval(1.42), sin(ComplexDouble(0,1)*1.42)))
	return 15;
  if(!DBL_EQUAL(Fcos.Eval(1.42), cos(1.42)))
	return 16;
  if(!DBL_EQUAL(Fcosp.Eval(1.42), cos(1.42)))
	return 17;
  if(!DBL_EQUAL(Fcosip.Eval(1.42), cos(ComplexDouble(0,1)*1.42)))
	return 18;
  if(!DBL_EQUAL(Fexp.Eval(1.42), exp(1.42)))
	return 19;
  if(!DBL_EQUAL(Fexpp.Eval(1.42), exp(1.42)))
	return 20;
  if(!DBL_EQUAL(Fexpm.Eval(1.42), exp(-1.42)))
	return 21;
  if(!DBL_EQUAL(Fexpip.Eval(1.42), exp(ComplexDouble(0,1)*1.42)))
	return 22;
  if(!DBL_EQUAL(Fexpim.Eval(1.42), exp(-ComplexDouble(0,1)*1.42)))
	return 23;

  return 0;
}

int BasisFunctionTest::runUnitTests() const
{
  cout << "Running unit tests on BasisFunction...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;

  cout << "done" << endl;
  return 0;
}
