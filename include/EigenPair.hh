#ifndef EigenPair_hh
#define EigenPair_hh 1

#include <vector>
#include <cstring>
#include "Globals.hpp"


using namespace std;


class EigenPair
{
public:
  ComplexDouble Eigenvalue;
  vector<ComplexDouble> Eigenvector;  
  bool operator<(const EigenPair &other) const;

  EigenPair();
  EigenPair(size_t size);
private:
  static int ComplexFirst(const ComplexDouble & c1, 
						   const ComplexDouble & c2
						   ); ///Returns 1 if c1 is smaller than c2 according to some arbitrary but well-defined ordering (sort by real part followed by imaginary part.). Returns -1 if c2 is smaller, and 0 if they are equal.
};
#endif
