#ifndef PrecomputedInteractionEvaluator_hh
#define PrecomputedInteractionEvaluator_hh 1


#include <vector>
#include <iostream>
#include "CompositeBasisFunction.hh"
#include "InteractionProperties.hh"
#include "LegendreRule.hh"
#include "UnorderedTuple.hh"
#include "HarmonicBasisFunction.hh"
#include "SpecificUnits.hh"
#include "HermiteRule.hh"
#include "HermiteEvaluator.hh"
#include "VerbosePrinter.hh"
#include "MultiTasker.hpp"
#include "ElementComputationWorkerData.hh"
using namespace std;

#define FILE_VALIDATION_NUMBER 1357955





/** Since it's expensive to call Eval() on a CompositeBasisFunction object,
	this object is an intermediate which has Eval()uated all the required points on the CompositeBasisFunctions already. A call to Eval() here instead evaluates the interaction element directly, which should be faster and also parallellalizable without issues.

 */
class PrecomputedInteractionEvaluator
{
public:
  PrecomputedInteractionEvaluator(vector<CompositeBasisFunction *> * myBasisFunctions, ///Basis functions to use.
								  const InteractionProperties * myInteractionProperties, ///Interaction properties.
								  const HarmonicBasisFunction * myHarmonicBasisFunction, ///Harmonic basis function.
								  const SpecificUnits * mySpecificUnits, ///specific units.
								  VerbosePrinter * myPrinter = NULL, ///only used in the constructor.
								  uint numberOfThreads = 10///Number of threads to use for init.
								  );

  PrecomputedInteractionEvaluator(const InteractionProperties * myInteractionProperties,
								  VerbosePrinter * myPrinter = NULL
								  ); ///Initialize from file.
  
  ~PrecomputedInteractionEvaluator();
  
  ComplexDouble GetElement(uint a, ///First basis function.
						   uint b, ///Second basis function.
						   uint c, ///First basis function.
						   uint d ///Second basis function.
						   ) const;

  ComplexDouble GetEnergy(uint basis,
						  uint index
						  ) const;

  uint GetNumberOfBases() const;

  uint GetBasisElements(uint basisIndex
						) const;

  void SetNewCouplingCoefficient(double value);

  void DiskSave(const InteractionProperties * myInteractionProperties
				); ///Save object to file to be able to read it later on.


protected:
  static ComplexDouble ComputeSingleElement(uint a, 
											uint b, 
											uint c, 
											uint d, 
											const vector<vector<vector<vector<double> > > > *Vnnnn, 
											const vector<vector<vector<ComplexDouble> > > * PsiB, 
											double nmax
											); ///Computes an interaction matrix element (minus the coupling coefficient).

  PrecomputedInteractionEvaluator() ///For debugging purposes.
  {};

  void InitializeEnergies(vector<CompositeBasisFunction*> * myBasisFunctions
						  );

  static void ValidateInput(vector<CompositeBasisFunction*> * myBasisFunctions, 
							const InteractionProperties * myInteractionProperties
							);

  static vector<vector<vector<ComplexDouble > > > * ComputePsiB(vector<CompositeBasisFunction*> * myBasisFunctions, 
																uint integralPrecision, 
																const HarmonicBasisFunction * myHarmonicBasisFunction, 
																double bFactor, 
																VerbosePrinter * myPrinter,
																uint nmax
																);

  static vector<vector<vector<vector<double> > > > * ComputeVnnnn(uint integralPrecision, 
																  const HarmonicBasisFunction * myHarmonicBasisFunction, 
																  double bFactor, 
																  VerbosePrinter * myPrinter,
																  uint nmax
																  );

  void ComputeElements(vector<vector<vector<ComplexDouble> > > * PsiB, 
					   vector<vector<vector<vector<double> > > > * Vnnnn, 
					   VerbosePrinter * myPrinter,
					   uint nmax,
					   uint numberOfThreads
					   );

  static bool ComputeElementsDoWork(ElementComputationWorkerData * w
									);

private:
  vector< vector <vector <vector <ComplexDouble> > > > elements; ///Contains the precomputed interaction elements. We should have symmetry between the indices 1<->3 and 2<->4.


  vector<vector<ComplexDouble> > Energies; ///Contains energies corresponding to the basis functions. Outer layer is basis function number, inner layer is energy for that index. Not in general sorted in any way.

  double couplingCoefficient; ///Coupling coefficient.
};

#endif
