#ifndef ExternalLauncher_hh
#define ExternalLauncher_hh 1

#include "ComputeConfig.hh"
#include <stdlib.h>
#include <unistd.h>


class ExternalLauncher
{
public:
  ExternalLauncher(ComputeConfig * _config);
  void Launch() const;

private:
  ComputeConfig * config;
};

#endif
