#ifndef HermiteEvaluator_hh
#define HermiteEvaluator_hh 1

#include <vector>
#include "RLMacros.hpp"
#include "RLException.hh"
#include <iostream>
#include <iomanip>
using namespace std;


class HermiteEvaluator
{
public:
  static void Init(uint _n,
				   bool _autoResize = false
				   );

  static void DeInit(); ///Free allocated memory.

  static double HermiteH(uint n,
						 double x
						 );

protected:
  static void FillTable(uint nMax
						);

  static void DumpTable(); ///Use for debugging purposes.

private:
  static bool autoResize; ///Tells if we should automatically expand the coefficient table or throw an exception. Pros of auto-resize: convenience. Cons: not thread safe, the function will have side-effects.
  static vector<vector<double> > polyCoefficients;
  HermiteEvaluator() {} ///Prevent instantiation.
};

#endif
