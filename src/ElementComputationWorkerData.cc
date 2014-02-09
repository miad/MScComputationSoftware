#include "ElementComputationWorkerData.hh"

ElementComputationWorkerData::ElementComputationWorkerData(uint _a, uint _b, uint _nmax, vector<vector<vector<vector<ComplexDouble> > > > * _elements, const vector<vector<vector<ComplexDouble> > > * _PsiB, const vector<vector<vector<vector<double> > > > * _Vnnnn)
  :a(_a), b(_b), nmax(_nmax), elements(_elements), PsiB(_PsiB), Vnnnn(_Vnnnn)
{ 
  
}

ElementComputationWorkerData::~ElementComputationWorkerData()
{ 
  
}
