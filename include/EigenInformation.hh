#ifndef EigenInformation_hh
#define EigenInformation_hh 1


#include <vector>
#include "Globals.hpp"

class EigenInformation
{
public:
  EigenInformation();
  ~EigenInformation();
  

  vector<vector<ComplexDouble> > Eigenvectors;
  vector<ComplexDouble> Eigenvalues;

private:

};
#endif
