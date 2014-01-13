#include "HermiteEvaluatorTest.hh"

int HermiteEvaluatorTest::TestCase1() const
{
  HermiteEvaluator::DeInit();
  HermiteEvaluator::Init(0, false);
  try
	{
	  double vask = HermiteEvaluator::HermiteH(3.3, 1);
	  (void)vask;
	  return 1;
	}
  catch(exception &ex)
	{

	}
  ///Evaluate in the arbitrary point 4.4444.
  for(uint i = 0; i<2; ++i)
	{
	  cout.precision(25);
	  HermiteEvaluator::Init(20,i);
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(0,4.4444), 1L))
		return 100*i+10;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(1,4.4444), 8.8888L))
		return 100*i+11;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(2,4.4444), 77.010765440000000000000000000L))
		return 100*i+12;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(3,4.4444), 648.97809184307200000000000000L))
		return 100*i+13;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(4,4.4444), 5306.5718701346983936000000000L))
		return 100*i+14;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(5,4.4444), 41977.231304508731081031680000L))
		return 100*i+15;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(6,4.4444), 320061.49491817022489707439718L))
		return 100*i+16;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(7,4.4444), 2.3412358403745267220927347417E6L))
		return 100*i+17;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(8,4.4444), 1.6329916209066709978778859011E7L))
		  return 100*i+18;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(9,4.4444), 1.0769358575315974410588576611E8L))
		return 100*i+19;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(10,4.4444), 6.6332825327948555379037793562E8L))
		return 100*i+20;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(11,4.4444), 3.7423204626874963084141960719E9L))
		return 100*i+21;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(12,4.4444), 1.8671516556587935002843791460E10L))
		return 100*i+22;
	  if(!DBL_EQUAL_HPREC(HermiteEvaluator::HermiteH(13,4.4444), 7.61516852636989252513371878056194625490940605422090156798640E10L))
		return 100*i+23;
	  if(!DBL_EQUAL_MPREC(HermiteEvaluator::HermiteH(14,4.4444), 1.9143766950068069670014741700E11L))
		return 100*i+24;
	  if(!DBL_EQUAL_MPREC(HermiteEvaluator::HermiteH(15,4.4444), -4.305960307259193302091708983E11L))
		return 100*i+25;
	  if(!DBL_EQUAL_MPREC(HermiteEvaluator::HermiteH(16,4.4444), -9.5706120829369726433677007909E12L))
		return 100*i+26;

	}
  return 0;
}

int HermiteEvaluatorTest::runUnitTests() const
{
  cout << "Running unit tests on HermiteEvaluator...";
  cout << flush;
  int code1 = TestCase1();
  if(code1)
    return code1;

  cout << "done" << endl;
  return 0;
}
