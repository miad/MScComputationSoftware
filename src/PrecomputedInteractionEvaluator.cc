#include "PrecomputedInteractionEvaluator.hh"


PrecomputedInteractionEvaluator::PrecomputedInteractionEvaluator(vector<CompositeBasisFunction*> * myBasisFunctions, const InteractionProperties * myInteractionProperties, const HarmonicBasisFunction * myHarmonicBasisFunction, const SpecificUnits * mySpecificUnits, VerbosePrinter * myPrinter, uint numberOfThreads)
{
  if(myBasisFunctions == NULL || myInteractionProperties == NULL)
	{
	  throw RLException("PrecomputedInteractionEvaluator: invalid basis functions or interaction properties object: was NULL.");
	}

  if(myBasisFunctions->size() != 2)
	{
	  throw RLException("PrecomputedInteractionEvaluator only implemented for 2 basis functions.");
	}

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



  vector<vector<vector<ComplexDouble > > > PsiB;
  ///Outermost: composite basis function number.
  ///middle: element in the composite basis function
  ///inner: Hermite polynomial number.


  double couplingCoefficient; ///Coupling coefficient.
  uint nmax = myInteractionProperties->GetNMax();

  HermiteEvaluator::Init(nmax+10, false);



  couplingCoefficient = myInteractionProperties->GetCouplingCoefficient();


  double bFactor = 	mySpecificUnits->GetMassOverLambda2() * myHarmonicBasisFunction->GetOmega() / mySpecificUnits->GetHbarTimesLambda();


  //GHpoints = HermiteRule::GetRule(precision, xmin, MASS*omega/HBAR);
  vector<pair<double, double> > myHalfHermiteRule = 
	HermiteRule::GetRule(myInteractionProperties->GetPrecision(),
						 myHarmonicBasisFunction->GetXmin(), bFactor / 2.0);

  ///Here we can probably speed up things by using an intermediate cache to build up the solution in two steps. If required.
  
  ulong ele = 2 * myBasisFunctions->at(0)->GetSize() * myBasisFunctions->at(1)->GetSize() * nmax;
  
  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Interaction: Precomputing H-Psi values: %ld values to compute.\n", ele);
	  myPrinter->Print(4, "Progress: 00%%");
	}
  
  PsiB.resize(myBasisFunctions->size());
  
  uint pSteps1 = myBasisFunctions->at(0)->GetSize() +  myBasisFunctions->at(1)->GetSize();
  uint step = 0;

  for(uint i = 0; i<myBasisFunctions->size(); ++i)
	{
	  PsiB[i].resize(myBasisFunctions->at(i)->GetSize());
	  for(uint j = 0; j<myBasisFunctions->at(i)->GetSize(); ++j)
		{
		  if(myPrinter != NULL)
			{
			  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
			  myPrinter->Print(4, "Progress: %02i%%", 100*(step++)/pSteps1);
			}
		  
		  PsiB[i][j].resize(nmax);
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
			  PsiB[i][j][n] = ComplexDouble(real(sum), imag(sum));
			}
		}
	}


  vector<pair<double, double> > myDoubleHermiteRule = 
	HermiteRule::GetRule(myInteractionProperties->GetPrecision(),
						 myHarmonicBasisFunction->GetXmin(), bFactor * 2.0);


  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	  myPrinter->Print(3, "Interaction: Precomputing nnnn-values. %ld values to compute.\n", nmax*(nmax-1)*(nmax-2)*(nmax-3));
	  myPrinter->Print(4, "Progress: 00%%");
	}



  vector<vector<vector<vector<double> > > >  Vnnnn(nmax, vector<vector<vector<double> > >(nmax, vector<vector<double> >(nmax, vector<double>(nmax, 0.0) ) ) );
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
				  Vnnnn[n1][n2][n3][n4] = sum;
				}
			}
		}
	}

  if(myPrinter != NULL)
    {
	    myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
    }

  elements.resize(myBasisFunctions->at(0)->GetSize(),
				  vector< vector< vector< ComplexDouble > > > (myBasisFunctions->at(1)->GetSize(),
															   vector<vector<ComplexDouble> >(myBasisFunctions->at(0)->GetSize(), 
															   vector<ComplexDouble>(myBasisFunctions->at(1)->GetSize(), 
																					 0.0))));

  ulong nCross = pow(myBasisFunctions->at(0)->GetSize() * myBasisFunctions->at(1)->GetSize(), 2);

  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Interaction: Precomputing cross-values. %ld values to compute.\n", nCross);
	  myPrinter->Print(4, "Progress: 00.00%%");
	}


  ///TODO: can be simplified using symmetry.
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

				elements[a][b][c][d] = couplingCoefficient * ComputeElement(a, b, c, d, Vnnnn, PsiB, nmax);
				///Symmetry, should save us some computing time.
				elements[c][b][a][d] = elements[a][b][c][d];
				elements[c][d][a][b] = elements[a][b][c][d];
				elements[a][d][c][b] = elements[a][b][c][d];
			  }
		  }
	  }
  }
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	}

}

uint PrecomputedInteractionEvaluator::Permutations(uint a, uint b, uint c, uint d)
{
  short k = (a==b) + (a==c) + (a==d) + (b==c) + (b==d) + (c==d);
  switch(k)
	{
	case 0:
	  return 24;
	case 1:
	  return 12;
	case 2:
	  return 6;
	case 3:
	  return 4;
	case 6:
	  return 1;
	default:
	  throw RLException("Code design error.");
	}
}

ComplexDouble PrecomputedInteractionEvaluator::ComputeElement(uint a, uint b, uint c, uint d, vector<vector<vector<vector<double> > > > &Vnnnn, vector<vector<vector<ComplexDouble> > > & PsiB, double nmax)
{
  ComplexDouble sum = 0.0;
  for(uint n1 = 0; n1 < nmax; ++n1)
	{
	  for(uint n2 = n1; n2 < nmax; ++n2)
		{
		  for(uint n3 = n2; n3 < nmax; ++n3)
			{
			  for(uint n4 = n3; n4 < nmax; ++n4)
				{
				  sum += 
					PsiB[0][a][n1] * 
					PsiB[1][b][n2] * 
					PsiB[0][c][n3] * 
					PsiB[1][d][n4]
					* Vnnnn[n1][n2][n3][n4] * 
					(double)Permutations(n1, n2, n3, n4);
				}
			}
		}
	}
  return sum;
}







uint PrecomputedInteractionEvaluator::GetNumberOfBases() const
{
  return 2; ///for now. Function for portability.
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
