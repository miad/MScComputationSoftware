#include "PrecomputedInteractionEvaluator.hh"


PrecomputedInteractionEvaluator::PrecomputedInteractionEvaluator(vector<CompositeBasisFunction*> * myBasisFunctions, const InteractionProperties * myInteractionProperties, HarmonicBasisFunction * myHarmonicBasisFunction)
{

  if(myBasisFunctions == NULL || myInteractionProperties == NULL)
	{
	  throw RLException("PrecomputedInteractionEvaluator: invalid basis functions or interaction properties object: was NULL.");
	}
















  myLegendreRule = 
	LegendreRule::GetRule(
						  myInteractionProperties->GetPrecision(), 
						  myInteractionProperties->GetLowerIntegrationLimit(), 
						  myInteractionProperties->GetUpperIntegrationLimit()
						  );
  
  couplingCoefficient = myInteractionProperties->GetCouplingCoefficient();
  
  functionValues.resize(myBasisFunctions->size());
  Energies.resize(myBasisFunctions->size());

  ///Precompute a cache.
  for(uint i = 0; i<myBasisFunctions->size(); ++i)
	{
	  functionValues[i].resize(myBasisFunctions->at(i)->GetNumberOfParameters());
	  Energies.at(i).resize(myBasisFunctions->at(i)->GetNumberOfParameters());
	  for(uint j = 0; j<myBasisFunctions->at(i)->GetNumberOfParameters(); ++j)
		{
		  Energies[i][j] = myBasisFunctions->at(i)->GetE(j);
		  functionValues[i][j].resize(myLegendreRule.size(), 0.0);
		  for(uint k = 0; k<myLegendreRule.size(); ++k)
			{
			  functionValues[i][j][k] = myBasisFunctions->at(i)->Eval(myLegendreRule[k].first, j);
			}
		}
	}

}

uint PrecomputedInteractionEvaluator::GetNumberOfBases() const
{
  return functionValues.size();
}

uint PrecomputedInteractionEvaluator::GetBasisElements(uint basisIndex) const
{
  return functionValues.at(basisIndex).size();
}

ComplexDouble PrecomputedInteractionEvaluator::GetEnergy(uint basis, uint index) const
{
  return Energies.at(basis).at(index);
}


PrecomputedInteractionEvaluator::~PrecomputedInteractionEvaluator()
{

}


ComplexDouble PrecomputedInteractionEvaluator::GetElement(uint a, uint b, uint c, uint d) const
{
  if( DBL_EQUAL(couplingCoefficient, 0.0) )
	return 0.0;

  ComplexDouble sum = 0.0;
  for(uint i = 0; i<myLegendreRule.size(); ++i)
	{
	  sum += myLegendreRule[i].second * 
		functionValues.at(0).at(a).at(i) *
		functionValues.at(1).at(b).at(i) *
		functionValues.at(0).at(c).at(i) *
		functionValues.at(1).at(d).at(i);
	}
  
  return sum * couplingCoefficient;
}
