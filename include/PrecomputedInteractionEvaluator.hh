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
using namespace std;







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


protected:
  static ComplexDouble ComputeElement(uint a, 
									  uint b, 
									  uint c, 
									  uint d, 
									  vector<vector<vector<vector<double> > > > &Vnnnn, 
									  vector<vector<vector<ComplexDouble> > > & PsiB, 
									  double nmax
									  ); ///Computes an interaction matrix element (minus the coupling coefficient).

  static uint Permutations(uint a, 
						   uint b, 
						   uint c, 
						   uint d
						   ); ///Count number of permutations.

  PrecomputedInteractionEvaluator() {};


private:

  vector< vector <vector <vector <ComplexDouble> > > >elements;
  vector<vector<ComplexDouble> > Energies; 
  ///Outermost: basis function
  ///middle: basis parameter ID
  ///innermost: legendre x-point corresponding

};
#endif
