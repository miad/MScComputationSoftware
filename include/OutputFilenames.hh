#ifndef OutputFilenames_hh
#define OutputFilenames_hh 1

#include <string>
#include <map>

#include "RLException.hh"

using namespace std;

/*
  Contains the output filenames. Essentially a map, but which will throw exceptions if attempting to access a non-existent name.
 */
class OutputFilenames
{
public:
  OutputFilenames(); ///Constructor, initializes the empty object.
  void Add(string key, 
		   string value
		   ); ///Adds a key/value pair.
  bool Contains(string key
		   ) const; /// Check if the key is here.
  string Get(string key ///The key in question.
			 ) const; ///Returns the value for a key. Throws an exception if there is no value for that key.
  void Clear(); ///Clears the container and sets it to zero contained elements.
private:
  map<string, string> filenames; ///Underlying map.
};


#endif
