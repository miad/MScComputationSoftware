#ifndef EigenInformation_hh
#define EigenInformation_hh 1


#include <vector>
#include <algorithm>
#include "Globals.hpp"
#include "EigenPair.hh"

class EigenInformation
{
public:
  EigenInformation();
  ~EigenInformation();

  vector<EigenPair> EigenPairs;
};
#endif
