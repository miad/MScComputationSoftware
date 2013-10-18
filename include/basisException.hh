#include <exception>
#include <string>

#ifndef basisException_hh
#define basisException_hh 1

using namespace std;

class basisException : public exception
{
  
private:
  string msg; ///  The reason of the exception to be thrown.

public:
  basisException(string m="Exception!" /// The description of the reason for the exception to be thrown.
		    ); /// Constructor.

  ~basisException() throw(); /// Destructor.
  const char* what() const throw(); /// Retrieves the reason of the exception.
};

#endif
