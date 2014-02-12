#include "PrecomputedInteractionEvaluator.hh"

PrecomputedInteractionEvaluator::PrecomputedInteractionEvaluator(vector<CompositeBasisFunction*> * myBasisFunctions, const InteractionProperties * myInteractionProperties, const HarmonicBasisFunction * myHarmonicBasisFunction, const SpecificUnits * mySpecificUnits, VerbosePrinter * myPrinter, uint numberOfThreads)
{

  ValidateInput(myBasisFunctions, myInteractionProperties);
  InitializeEnergies(myBasisFunctions);

  
  uint nmax = myInteractionProperties->GetNMax();
  vector<double> bFactor;
  for(uint i = 0; i<2; ++i)
	bFactor.push_back(mySpecificUnits->GetMassOverLambda2() * myHarmonicBasisFunction->GetOmega(i) / mySpecificUnits->GetHbarTimesLambda());

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

  
  ComputeElements(PsiB, Vnnnn, myPrinter, nmax, numberOfThreads);

  couplingCoefficient = myInteractionProperties->GetCouplingCoefficient();
  
  delete PsiB; PsiB = NULL;
  delete Vnnnn; Vnnnn = NULL;
  if(myPrinter!=NULL)
	myPrinter->Print(3, "Optionally saving to disk...");
  DiskSave(myInteractionProperties);
  if(myPrinter != NULL)
	myPrinter->Print(3,"done\n");
}

void PrecomputedInteractionEvaluator::SetNewCouplingCoefficient(double value)
{
  couplingCoefficient = value;
}

ComplexDouble PrecomputedInteractionEvaluator::ComputeSingleElement(uint a, uint b, uint c, uint d, const vector<vector<vector<vector<double> > > > * Vnnnn, const vector<vector<vector<ComplexDouble> > > * PsiB, double nmax)
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
  return elements.at(a).at(b).at(c).at(d) * couplingCoefficient;
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

vector<vector<vector<ComplexDouble > > > * PrecomputedInteractionEvaluator::ComputePsiB(vector<CompositeBasisFunction*> * myBasisFunctions, uint integralPrecision, const HarmonicBasisFunction * myHarmonicBasisFunction, vector<double> & bFactor, VerbosePrinter * myPrinter, uint nmax)
{

  ///Outermost: composite basis function number.
  ///middle: element in the composite basis function
  ///inner: Hermite polynomial number.
  vector<vector<vector<ComplexDouble > > > * PsiB = new vector<vector<vector<ComplexDouble> > >();

  vector<vector<pair<double, double> > > myHalfHermiteRule;
  for(uint i = 0; i<2; ++i)
	myHalfHermiteRule.push_back(
								HermiteRule::GetRule(integralPrecision, myHarmonicBasisFunction->GetXmin(i), bFactor.at(i) / 2.0)
								);


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
			  for(vector<pair<double, double> > ::const_iterator it = myHalfHermiteRule.at(i).begin(); it!=myHalfHermiteRule.at(i).end(); ++it)
				{
				  ComplexDouble s = myBasisFunctions->at(i)->Eval(it->first, j);
				  complex<long double> rs = complex<long double>(real(s), imag(s));
				  sum += (complex<long double>)it->second * 
					myHarmonicBasisFunction->EvalNonExponentPart(n, it->first, i) * rs;
				}
			  PsiB->at(i).at(j).at(n) = ComplexDouble(real(sum), imag(sum));
			}
		}
	}
  return PsiB;
}



vector<vector<vector<vector<double> > > > * PrecomputedInteractionEvaluator::ComputeVnnnn(uint integralPrecision, const HarmonicBasisFunction * myHarmonicBasisFunction, vector<double> & bFactor, VerbosePrinter * myPrinter, uint nmax)
{


  double x1 = myHarmonicBasisFunction->GetXmin(0);
  double x2 = myHarmonicBasisFunction->GetXmin(1);
  double b1 = bFactor.at(0);
  double b2 = bFactor.at(1);

  double bTilde = (b1+b2);
  double xTilde = (b1*x1+b2*x2)/(b1+b2);
  double preFExp = (b1*x1*x1+b2*x2*x2) - xTilde*xTilde*bTilde;
  double integralPrefactor = exp(-preFExp);



  vector<pair<double, double> > myDoubleHermiteRule = HermiteRule::GetRule(integralPrecision, 
											 xTilde, bTilde);


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
					  sum += integralPrefactor * it->second * 
						myHarmonicBasisFunction->EvalNonExponentPart(n1, it->first, 0) *
						myHarmonicBasisFunction->EvalNonExponentPart(n2, it->first, 1) *
						myHarmonicBasisFunction->EvalNonExponentPart(n3, it->first, 0) *
						myHarmonicBasisFunction->EvalNonExponentPart(n4, it->first, 1);
						
					}
				  uint vals[] = {n1, n2, n3, n4};
				  do
					{
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


void PrecomputedInteractionEvaluator::ComputeElements(vector<vector<vector<ComplexDouble> > > * PsiB, vector<vector<vector<vector<double> > > > * Vnnnn, VerbosePrinter * myPrinter, uint nmax, uint numberOfThreads)
{

  elements.resize(PsiB->at(0).size(), 
				  vector< vector< vector< ComplexDouble > > > (PsiB->at(1).size(),
				  vector<vector<ComplexDouble> >(PsiB->at(0).size(), 
				  vector<ComplexDouble>(PsiB->at(1).size(), 
				  0.0))));
  
  ulong nCross = pow(PsiB->at(0).size() * PsiB->at(1).size(), 2);

  MultiTasker<ElementComputationWorkerData*, bool> * myMultiTasker = new MultiTasker<ElementComputationWorkerData*, bool>(ComputeElementsDoWork, numberOfThreads);
  myMultiTasker->RegisterListener(myPrinter);

  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Queuing up cross-value work.\n");
	}

  ulong workElements = 0;
  for(uint a = 0; a<elements.size(); ++a)
  {
	for(uint b = 0; b<elements.at(a).size(); ++b)
	  {
		++workElements;
		myMultiTasker->AddInput(new ElementComputationWorkerData(a, b, nmax, &elements, PsiB, Vnnnn));
	  }
  }

  if(myPrinter != NULL)
	{
	  myPrinter->Print(3, "Interaction: Precomputing cross-values. %ld values to compute.\n", nCross);
	  myPrinter->Print(4, "Progress: 00.00%%");
	}
  
  
  myMultiTasker->LaunchThreads();
  ulong progress = 0;
  while(myMultiTasker->GetInputMinusOutput())
	{
	  myMultiTasker->GetOutput();
	  if(myPrinter != NULL)
		{
		  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
		  myPrinter->Print(4, "Progress: %04.02f%%", (double)100*((double)progress++/((double)workElements)));
		}
	}
  
  myMultiTasker->PauseUntilOutputIsGenerated();
  myMultiTasker->DestroyThreads();
  
  
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	}
  delete myMultiTasker; 
  myMultiTasker = NULL;
}



bool PrecomputedInteractionEvaluator::ComputeElementsDoWork(ElementComputationWorkerData * w)
{
  uint a = w->a;
  uint b = w->b;
  uint nmax = w->nmax;
  for(uint c = 0; c<=a; ++c)
	{
	  for(uint d = 0; d<=b; ++d)
		{
		  ComplexDouble val = PrecomputedInteractionEvaluator::ComputeSingleElement(a, b, c, d, w->Vnnnn, w->PsiB, nmax);
		  uint elPtr1[] = {a, c};
		  uint elPtr2[] = {b, d};
		  do
			{
			  do
				{
				  w->elements->at(elPtr1[0])[elPtr2[0]][elPtr1[1]][elPtr2[1]] = val;
				}
			  while(prev_permutation(elPtr2, elPtr2+2));
			}
		  while(prev_permutation(elPtr1, elPtr1+2));
		}
	}
  delete w;
  return true; ///We need to return _something_ in order to count the number of returned values for progress bar.
}


PrecomputedInteractionEvaluator::PrecomputedInteractionEvaluator(const InteractionProperties * myInteractionProperties, VerbosePrinter * myPrinter)
{
  couplingCoefficient = myInteractionProperties->GetCouplingCoefficient();
  string filename = myInteractionProperties->GetCacheFile();
  if(filename.empty()) ///Don't save if filename is empty.
	{
	  throw RLException("Cannot read from non-existing input file.");
	}
  
  FILE * fin = fopen(filename.c_str(), "rb");
  if(fin == NULL)
	{
	  throw RLException("Could not open interaction output file.");
	}
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "Using file '%s'.\n", filename.c_str());
	}
  ulong precision;
  uint nmax;
  uint N[2];
  if(fread(&precision, sizeof(precision), 1, fin) != 1)
	throw RLException("File read of precision failed.");
  if(fread(&nmax, sizeof(nmax), 1, fin) != 1)
	throw RLException("Failed read of nmax");

  if(myInteractionProperties->GetNMax() != nmax || myInteractionProperties->GetPrecision() != precision)
	{
	  throw RLException("Interaction data in cache file does not match some config settings.");
	}


  if(fread(&N[0], sizeof(N[0]), 2, fin) != 2)
	throw RLException("Failed read of N(1 or 2)");
  Energies.resize(2); 
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "Loading energies...", filename.c_str());
	}

  for(uint i = 0; i<2; ++i)
	{
	  Energies[i].resize(N[i]);
	  if(fread(&Energies[i][0], sizeof(Energies[i][0]), N[i], fin) != N[i])
		throw RLException("Failed read.");
	}
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "done\n", filename.c_str());
	  myPrinter->Print(4, "Mapping RAM memory...", filename.c_str());
	}




  elements.resize(N[0], 
				  vector< vector< vector< ComplexDouble > > > (N[1],
				  vector<vector<ComplexDouble> >(N[0], 
				  vector<ComplexDouble>(N[1], 
				  0.0))));

  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "done\n");
	  myPrinter->Print(3, "Loading elements from file.\n");
	  myPrinter->Print(4, "Progress: 00.00%%");
	}


  for(uint a = 0; a<N[0]; ++a)
	{
	  for(uint b = 0; b<N[1]; ++b)
		{
		  if(myPrinter != NULL)
			{
			  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
			  myPrinter->Print(4, "Progress: %04.02f%%", (double)100*((double)(a*N[1]+b)/((double)N[0]*N[1])));
			}

		  for(uint c = 0; c<N[0]; ++c)
			{
			  if(fread(&elements[a][b][c][0], sizeof(elements[a][b][c][0]), N[1], fin) != N[1])
				throw RLException("File read failed : a=%d, b=%d, c=%d, dsize=%d", a, b, c, N[1]);
			}
		}
	}

  ulong validateNo = 0;
  if(fread(&validateNo, sizeof(validateNo), 1, fin) != 1)
	throw RLException("File read of validateNo failed.");
  if(validateNo != FILE_VALIDATION_NUMBER)
	{
	  throw RLException("File consistency validation failed.");
	}

  fclose(fin); fin = NULL;
  if(myPrinter != NULL)
	{
	  myPrinter->Print(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	}

}

void PrecomputedInteractionEvaluator::DiskSave(const InteractionProperties * myInteractionProperties)
{
  string filename = myInteractionProperties->GetCacheFile();
  if(filename.empty()) ///Don't save if filename is empty.
	return; 

  FILE * fout = fopen(filename.c_str(), "wb");
  if(fout == NULL)
	{
	  throw RLException("Could not open interaction output file.");
	}
  ///Save interaction properties for data validation on readback.
  ulong precision = myInteractionProperties->GetPrecision();
  uint nmax = myInteractionProperties->GetNMax();
  if(fwrite(&precision, sizeof(precision), 1, fout) != 1)
	throw RLException("File write failed.");
  if(fwrite(&nmax, sizeof(nmax), 1, fout) != 1)
	throw RLException("File write failed.");

  uint N[2];
  for(uint i = 0; i<2; ++i)
	N[i] = Energies[i].size();

  ///Save the two primary numbers: the number of basis elements for each particle.
  if(fwrite(&N[0], sizeof(N[0]), 2, fout) != 2)
	throw RLException("File write failed.");
  
  ///Save individual particle energies.
  for(uint i = 0; i<2; ++i)
	{
	  if(fwrite(&Energies[i][0], sizeof(Energies[i][0]), N[i], fout) != Energies[i].size())
		throw RLException("File write failed.");
	  
	}
  
  if(elements.size() != N[0])
	throw RLException("Size mismatch.");
  for(uint a = 0; a<N[0]; ++a)
	{
	  if(elements[a].size() != N[1])
		{
		  throw RLException("Size mismatch");
		}
	  for(uint b = 0; b<N[1]; ++b)
		{
		  if(elements[a][b].size() != N[0])
			{
			  throw RLException("Size mismatch");
			}
		  
		  for(uint c = 0; c<N[0]; ++c)
			{
			  if(elements[a][b][c].size() != N[1])
				{
				  throw RLException("Size mismatch");
				}
			  if(fwrite(&elements[a][b][c][0], sizeof(elements[a][b][c][0]), N[1], fout) != N[1])
				throw RLException("File write failed.");
			}
		}
	}

  ulong validateNo = FILE_VALIDATION_NUMBER;
  if(fwrite(&validateNo, sizeof(validateNo), 1, fout) != 1)
	throw RLException("File write failed.");
  
  fclose(fout); fout = NULL;
}
