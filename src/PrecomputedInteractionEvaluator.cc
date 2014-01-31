#include "PrecomputedInteractionEvaluator.hh"

PrecomputedInteractionEvaluator::PrecomputedInteractionEvaluator(vector<CompositeBasisFunction*> * myBasisFunctions, const InteractionProperties * myInteractionProperties, const HarmonicBasisFunction * myHarmonicBasisFunction, const SpecificUnits * mySpecificUnits, VerbosePrinter * myPrinter, uint numberOfThreads)
{

  ValidateInput(myBasisFunctions, myInteractionProperties);
  InitializeEnergies(myBasisFunctions);

  
  uint nmax = myInteractionProperties->GetNMax();
  double bFactor = 	mySpecificUnits->GetMassOverLambda2() * myHarmonicBasisFunction->GetOmega() / mySpecificUnits->GetHbarTimesLambda();
  HermiteEvaluator::Init(myInteractionProperties->GetNMax(), false);



  vector<vector<vector<ComplexDouble > > > * PsiB = ComputePsiB(myBasisFunctions, 
																myInteractionProperties->GetPrecision(), 
																myHarmonicBasisFunction,
																bFactor,
																myPrinter,
																nmax
																);

  vector<vector<vector<vector<double> > > > * Vnnnn = ComputeVnnnn(myInteractionProperties->GetPrecision(), 
																   myHarmonicBasisFunction,
																   bFactor,
																   myPrinter,
																   nmax
																   );

  ComputeElements(PsiB, Vnnnn, myPrinter, myInteractionProperties->GetCouplingCoefficient(), nmax);
  
  delete PsiB; PsiB = NULL;
  delete Vnnnn; Vnnnn = NULL;
}

ComplexDouble PrecomputedInteractionEvaluator::ComputeSingleElement(uint a, uint b, uint c, uint d, vector<vector<vector<vector<double> > > > * Vnnnn, vector<vector<vector<ComplexDouble> > > * PsiB, double nmax)
{
  ComplexDouble sum = 0.0;
  for(uint n1 = 0; n1 < nmax; ++n1)
	{
	  ComplexDouble ProdA = PsiB->at(0)[a][n1];
	  for(uint n2 = 0; n2 < nmax; ++n2)
		{
		  ComplexDouble ProdB = ProdA * PsiB->at(1)[b][n2];
		  for(uint n3 = 0; n3 < nmax; ++n3)
			{
			  ComplexDouble ProdC = ProdB * PsiB->at(0)[c][n3];
			  for(uint n4 = 0; n4 < nmax; ++n4)
				{
				  ComplexDouble ProdD = ProdC * PsiB->at(1)[d][n4];


				  sum += ProdD * Vnnnn->at(n1)[n2][n3][n4];
				}
			}
		}
	}
  return sum;
}



uint PrecomputedInteractionEvaluator::GetNumberOfBases() const
{
  return 2; ///for now. Function for portability in the future.
}

uint PrecomputedInteractionEvaluator::GetBasisElements(uint basisIndex) const
{
  switch(basisIndex)
	{
	case 0:
	  {
		return elements.size();
	  }
	case 1:
	  {
		return elements.at(0).size();
	  }
	default:
	  {
		throw RLException("Invalid basis index.");
	  }
	}
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
  return elements.at(a).at(b).at(c).at(d);
}



void PrecomputedInteractionEvaluator::InitializeEnergies(vector<CompositeBasisFunction*> * myBasisFunctions)
{
  Energies.resize(myBasisFunctions->size());
  
  ///Precompute a cache.
  for(uint i = 0; i<myBasisFunctions->size(); ++i)
	{
	  Energies.at(i).resize(myBasisFunctions->at(i)->GetSize());
	  for(uint j = 0; j<myBasisFunctions->at(i)->GetSize(); ++j)
		{
		  Energies[i][j] = myBasisFunctions->at(i)->GetE(j);
		}
	}
}

void PrecomputedInteractionEvaluator::ValidateInput(vector<CompositeBasisFunction*> * myBasisFunctions, const InteractionProperties * myInteractionProperties)
{
  if(myBasisFunctions == NULL || myInteractionProperties == NULL)
	{
	  throw RLException("PrecomputedInteractionEvaluator: invalid basis functions or interaction properties object: was NULL.");
	}
  
  if(myBasisFunctions->size() != 2)
	{
	  throw RLException("PrecomputedInteractionEvaluator only implemented for 2 basis functions.");
	}

}

vector<vector<vector<ComplexDouble > > > * PrecomputedInteractionEvaluator::ComputePsiB(vector<CompositeBasisFunction*> * myBasisFunctions, uint integralPrecision, const HarmonicBasisFunction * myHarmonicBasisFunction, double bFactor, VerbosePrinter * myPrinter, uint nmax)
{

  ///Outermost: composite basis function number.
  ///middle: element in the composite basis function
  ///inner: Hermite polynomial number.
  vector<vector<vector<ComplexDouble > > > * PsiB = new vector<vector<vector<ComplexDouble> > >();

  vector<pair<double, double> > myHalfHermiteRule = 
	HermiteRule::GetRule(integralPrecision,
						 myHarmonicBasisFunction->GetXmin(), bFactor / 2.0);


  ulong ele = 2 * myBasisFunctions->at(0)->GetSize() * myBasisFunctions->at(1)->GetSize() * nmax;
  
  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Interaction: Precomputing H-Psi values: %ld values to compute.\n", ele);
	  myPrinter->Print(4, "Progress: 00%%");
	}
  
  PsiB->resize(myBasisFunctions->size());
  
  uint pSteps1 = myBasisFunctions->at(0)->GetSize() +  myBasisFunctions->at(1)->GetSize();
  uint step = 0;

  for(uint i = 0; i<myBasisFunctions->size(); ++i)
	{
	  PsiB->at(i).resize(myBasisFunctions->at(i)->GetSize());
	  for(uint j = 0; j<myBasisFunctions->at(i)->GetSize(); ++j)
		{
		  if(myPrinter != NULL)
			{
			  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
			  myPrinter->Print(4, "Progress: %02i%%", 100*(step++)/pSteps1);
			}
		  
		  PsiB->at(i).at(j).resize(nmax);
		  for(uint n = 0; n<nmax; ++n)
			{
			  ///Summation should be done using extended precision.
			  complex<long double> sum = 0.0;
			  for(vector<pair<double, double> > ::const_iterator it = myHalfHermiteRule.begin(); it!=myHalfHermiteRule.end(); ++it)
				{
				  ComplexDouble s = myBasisFunctions->at(i)->Eval(it->first, j);
				  complex<long double> rs = complex<long double>(real(s), imag(s));
				  sum += (complex<long double>)it->second * 
					myHarmonicBasisFunction->EvalNonExponentPart(n, it->first) * rs;
				}
			  PsiB->at(i).at(j).at(n) = ComplexDouble(real(sum), imag(sum));
			}
		}
	}
  return PsiB;
}



vector<vector<vector<vector<double> > > > * PrecomputedInteractionEvaluator::ComputeVnnnn(uint integralPrecision, const HarmonicBasisFunction * myHarmonicBasisFunction, double bFactor, VerbosePrinter * myPrinter, uint nmax)
{
  
  vector<pair<double, double> > myDoubleHermiteRule = 
	HermiteRule::GetRule(integralPrecision, 
						 myHarmonicBasisFunction->GetXmin(), bFactor * 2.0);


  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	  myPrinter->Print(3, "Interaction: Precomputing nnnn-values. %ld values to compute.\n", nmax*(nmax-1)*(nmax-2)*(nmax-3));
	  myPrinter->Print(4, "Progress: 00%%");
	}

  vector<vector<vector<vector<double> > > >  * Vnnnn = new vector<vector<vector<vector<double> > > > (nmax, vector<vector<vector<double> > >(nmax, vector<vector<double> >(nmax, vector<double>(nmax, 0.0) ) ) );
  ///from outer to inner vectors: n1, n2, n3, n4.
  ///value: integral Psi_n1 Psi_n2 Psi_n1' Psi_n2' \dd x

  for(uint n1 = 0; n1 < nmax; ++n1)
	{
	  for(uint n2 = n1; n2 < nmax; ++n2)
		{
		  if(myPrinter != NULL)
			{
			  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
			  myPrinter->Print(4, "Progress: %02i%%", 100*(nmax*n1 + n2)/(nmax*nmax));
			}
		  for(uint n3 = n2; n3 < nmax; ++n3)
			{
			  for(uint n4 = n3; n4 < nmax; ++n4)
				{
				  long double sum = 0.0;
				  for(vector<pair<double, double> > ::const_iterator it = myDoubleHermiteRule.begin(); it!=myDoubleHermiteRule.end(); ++it)
					{
					  sum += it->second * 
						myHarmonicBasisFunction->EvalNonExponentPart(n1, it->first) *
						myHarmonicBasisFunction->EvalNonExponentPart(n2, it->first) *
						myHarmonicBasisFunction->EvalNonExponentPart(n3, it->first) *
						myHarmonicBasisFunction->EvalNonExponentPart(n4, it->first);
						
					}
				  uint vals[] = {n1, n2, n3, n4};
				  uint nPr = 0;
				  do
					{
					  ++nPr;
					  Vnnnn->at(vals[0]).at(vals[1]).at(vals[2]).at(vals[3]) = sum;
					}while(next_permutation(vals, vals+4));
				}
			}
		}
	}

  if(myPrinter != NULL)
    {
	    myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
    }
  return Vnnnn;
}


void PrecomputedInteractionEvaluator::ComputeElements(vector<vector<vector<ComplexDouble> > > * PsiB, vector<vector<vector<vector<double> > > > * Vnnnn, VerbosePrinter * myPrinter, double couplingCoefficient, uint nmax)
{

  elements.resize(PsiB->at(0).size(), 
				  vector< vector< vector< ComplexDouble > > > (PsiB->at(1).size(),
				  vector<vector<ComplexDouble> >(PsiB->at(0).size(), 
				  vector<ComplexDouble>(PsiB->at(1).size(), 
				  0.0))));
  
  ulong nCross = pow(PsiB->at(0).size() * PsiB->at(1).size(), 2);

  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Interaction: Precomputing cross-values. %ld values to compute.\n", nCross);
	  myPrinter->Print(4, "Progress: 00.00%%");
	}


  for(uint a = 0; a<elements.size(); ++a)
  {
	for(uint b = 0; b<elements.at(a).size(); ++b)
	  {
		if(myPrinter != NULL)
		  {
			myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
			myPrinter->Print(4, "Progress: %04.02f%%", (double)100*(a * elements.at(a).size() + b ) / (double)(elements.at(a).size() * elements.size() ));
		  }
		for(uint c = 0; c<=a; ++c)
		  {
			for(uint d = 0; d<=b; ++d)
			  {
				ComplexDouble val = couplingCoefficient * ComputeSingleElement(a, b, c, d, Vnnnn, PsiB, nmax);
				uint elPtr1[] = {a, c};
				uint elPtr2[] = {b, d};
				
				uint nPr = 0;
				do
				  {
					do
					  {
						elements[elPtr1[0]][elPtr2[0]][elPtr1[1]][elPtr2[1]] = val;
						++nPr;
					  }
					while(prev_permutation(elPtr2, elPtr2+2));
				  }
				while(prev_permutation(elPtr1, elPtr1+2));
				
			  }
		  }
	  }
  }

  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	}
}
