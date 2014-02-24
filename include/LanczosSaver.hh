#ifndef LanczosSaver_hh
#define LanczosSaver_hh 1
#include "RLMacros.hpp"
#include <string>
#include <cstring>
#include "RLException.hh"
#include "Globals.hpp"

#include <iostream>

#define FILE_VALIDATION_NUMBER_LCS 38198464



class LanczosSaver
{
public:
  static bool Save(const char * filename,
				   const CMatrix * matrix,
				   int * r
				   );
								 
private:
  LanczosSaver() { };
};
#endif
