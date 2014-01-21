#include "UnorderedTuple.hh"

UnorderedTuple::UnorderedTuple(vector<uint> _numbers)
  :numbers(_numbers)
{
  sort(numbers.begin(), numbers.end());
}

UnorderedTuple::UnorderedTuple(uint a, uint b, uint c, uint d)
{
  numbers.push_back(a); numbers.push_back(b); numbers.push_back(c); numbers.push_back(d);
  sort(numbers.begin(), numbers.end());
}

bool UnorderedTuple::operator==(const UnorderedTuple & rhs) const
{
  if(numbers.size() != rhs.numbers.size())
	return false;
  for(uint i = 0; i<numbers.size(); ++i)
	{
	  if(numbers[i] != rhs.numbers[i])
		return false;
	}
  return true;
}
