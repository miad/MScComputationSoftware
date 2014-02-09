#ifndef ElementComputationWorkerData_hh
#define ElementComputationWorkerData_hh 1


#include <vector>
#include "RLMacros.hpp"
#include "Globals.hpp"

class ElementComputationWorkerData
{
public:
  ElementComputationWorkerData(uint _a, uint _b, uint _nmax, vector<vector<vector<vector<ComplexDouble> > > > * _elments, const vector<vector<vector<ComplexDouble> > > * _PsiB, const vector<vector<vector<vector<double> > > > * _Vnnnn);
  ~ElementComputationWorkerData();

public:
  uint a;
  uint b;
  uint nmax;
  vector< vector <vector <vector <ComplexDouble> > > > * elements;
  const vector<vector<vector<ComplexDouble> > > * PsiB;
  const vector<vector<vector<vector<double> > > > * Vnnnn;

};
#endif
