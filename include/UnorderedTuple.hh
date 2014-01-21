#ifndef UnorderedTuple_hh
#define UnorderedTuple_hh 1
#include <vector>
#include "RLMacros.hpp"
#include <algorithm>
using namespace std;


class UnorderedTuple
{
public:
  UnorderedTuple(vector<uint> _numbers);
  UnorderedTuple(uint a, uint b, uint c, uint d);
  bool operator == (const UnorderedTuple & rhs) const;
private:
  vector<uint> numbers;
};
#endif
