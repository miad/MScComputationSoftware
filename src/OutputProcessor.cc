#include "OutputProcessor.hh"

OutputProcessor::OutputProcessor(ComputeConfig * _config)
  :config(_config), eigenData(NULL), myCompositeBasisFunctions(NULL)
{
  if(config == NULL)
	{
	  throw RLException("Cannot send a NULL config to OutputProcessor. ");
	}
}

OutputProcessor::~OutputProcessor()
{

}

void OutputProcessor::SetEigenInformation(EigenInformation * data)
{
  eigenData = data;
}

void OutputProcessor::SetCompositeBasisFunctions(vector<CompositeBasisFunction*> * value)
{
  myCompositeBasisFunctions = value;
}


void OutputProcessor::WritePostOutput() const
{
  if(eigenData == NULL)
	throw RLException("OutputProcessor error: eigenData was null.");

  WritePotentialToFile();
  WritePotentialPrecisionToFile();
  WriteKCurveToFile();
  WriteKFoundToFile();
  WriteEnergiesToFile();
  SavePrimaryEigenvectors();

  if(config->GetNumberOfParticles() == 1)
	{
	  WriteInterestingOneParticleWavefunctionsToFile();

	  WriteInterestingKPointsVerbosely(); ///Auto-ID:ed by filter algorithm
	  WriteInterestingKPointsToFile(); ///Auto-ID:ed by filter algorithm.
	}
  else if(config->GetNumberOfParticles() == 2)
	{
	  WriteInterestingRelativeTwoParticleWavefunctionsToFile();
	  WriteProductTwoParticleWavefunctionToFile();
	}
}

void OutputProcessor::SavePrimaryEigenvectors() const
{
  string fileName = config->GetOutputFilenames()->Get("PrimaryEigenvectors");
  if(fileName.empty() )
	{
	  vPrint(4, "Empty PrimaryEigenvectors filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving primary eigenvectors...");

  FILE * fout = AssuredFopen(fileName);

  fprintf(fout, "#Format: Re(eigenval) Im(Eigenval) Re(eigvect[0]) Im(eigvect[0]) etc.\n");

  for(ulong i = 0; i<eigenData->EigenPairs.size(); ++i)
	{
	  fprintf(fout, "%+13.10e %+13.10e", real(eigenData->EigenPairs[i].Eigenvalue), imag(eigenData->EigenPairs[i].Eigenvalue));
	  for(ulong j = 0; j<eigenData->EigenPairs[i].Eigenvector.size(); ++j)
		{
		  fprintf(fout, " %+13.10e %+13.10e", real(eigenData->EigenPairs[i].Eigenvector[j]), imag(eigenData->EigenPairs[i].Eigenvector[j]));
		}
	  fprintf(fout, "\n");
	}

  fclose(fout);
  fout = NULL; 

  vPrint(4, "done\n");
}


void OutputProcessor::WritePotentialToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("PotentialFile");
  if(fileName.empty() )
	{
	  vPrint(4, "Empty potential filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving potential...");

  FILE * fout = AssuredFopen(fileName);

  Potential * potential = config->GetPotential();

  double potentialLength = potential->GetMaxX() - potential->GetMinX();
  if(potentialLength < 0)
	throw RLException("The potential length was calculated to be < 0. Something is wrong.");
  
  double scaleFactor = 0.2; ///How much zero spacing to add at the end of the potential.

  vector<pair<double, double> >  plottingPoints = potential->GetPlottingPoints();

  fprintf(fout, "%+13.10e %+13.10e\n", potential->GetMinX() - scaleFactor*potentialLength, 0.);
  for(vector<pair<double, double> >::const_iterator it = plottingPoints.begin(); it!=plottingPoints.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", it->first, it->second);
	}
  fprintf(fout, "%+13.10e %+13.10e\n", potential->GetMaxX() + scaleFactor*potentialLength,0.0);
 
  fclose(fout);
  fout = NULL; 

  vPrint(4, "done\n");
}

void OutputProcessor::WritePotentialPrecisionToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("PotentialPrecisionFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty potential precision filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving potential precision points...");

  Potential * potential = config->GetPotential();

  FILE * fout = AssuredFopen(fileName);

  vector<pair<double, double> >  plottingPoints = potential->GetPrecisionPoints();
  for(vector<pair<double, double> >::const_iterator it = plottingPoints.begin(); it!=plottingPoints.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", it->first, it->second);
	}
  
  fclose(fout);
  fout = NULL; 

  vPrint(4, "done\n");
}

void OutputProcessor::WriteKCurveToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("KCurveFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty KCurve filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving K-curve...");
  ParametrizedCurve * toPrint = config->GetKCurve();
  FILE * fout = AssuredFopen(fileName);

  for(ulong i = 0; i<toPrint->GetNumberOfSegments(); ++i)
	{
	  const vector<pair<ComplexDouble, ComplexDouble> > * rule = toPrint->GetSegmentRule(i);
	  for(vector<pair<ComplexDouble, ComplexDouble> >::const_iterator it = rule->begin(); it != rule->end(); ++it)
		{
		  fprintf(fout, "%+13.10e %+13.10e\n", real(it->first), imag(it->first));
		}
	}
  
  fclose(fout);
  fout = NULL;
  
  vPrint(4, "done\n");
}


FILE * OutputProcessor::AssuredFopen(const string filename)
{
  FILE * fout = fopen(filename.c_str(), "w");
  if(fout == NULL)
	throw RLException("Could not open file '%s' for output.", filename.c_str());
  return fout;
}

void OutputProcessor::WriteKFoundToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("KFoundFile");

  if( fileName.empty() )
	{
	  vPrint(4, "Empty KFound filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving k-values for found eigenvalues...");

  FILE * fout = AssuredFopen(fileName);

  for(vector<EigenPair>::const_iterator it = eigenData->EigenPairs.begin(); it!=eigenData->EigenPairs.end(); ++it)
	{
	  ComplexDouble kToPrint = config->GetSpecificUnits()->EnergyToKValue(it->Eigenvalue);

	  ///If numerical stability is mean to us, then rotate.
	  if( (abs(imag(kToPrint)) > 1E2*abs(real(kToPrint))  && imag(kToPrint) < 0))
		{
		  kToPrint *= -1.0;
		}
	  fprintf(fout, "%+13.10e %+13.10e\n", real(kToPrint), imag(kToPrint));
	}
  
  fclose(fout);
  fout = NULL;
  vPrint(4, "done\n");
}


void OutputProcessor::WriteEnergiesToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("EnergyFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty Energy filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving found eigenvalues (energies)...");

  FILE * fout = AssuredFopen(fileName);


  for(vector<EigenPair>::const_iterator it = eigenData->EigenPairs.begin(); it!=eigenData->EigenPairs.end(); ++it)
	{
	  fprintf(fout, "%+13.10e %+13.10e\n", real(it->Eigenvalue), imag(it->Eigenvalue));
	}
  
  fclose(fout);
  fout = NULL;
  vPrint(4, "done\n");
}



void OutputProcessor::WriteInterestingKPointsVerbosely(bool forceFilter, ulong pid) const
{
  vector<ComplexDouble> printVector = FindInterestingKPoints(forceFilter);
  vector<ulong> idxVector = FindInterestingKPointIndex(forceFilter);

  for(ulong id = 0; id<printVector.size(); ++id)
	{
	  double RVu = real(config->GetSpecificUnits()->KValueToEnergy(printVector[id]) - config->GetPotential(pid)->Evaluate(config->GetHarmonicBasisFunction()->GetXmin(pid)) ); ///Energy in custom unit.

	  double RVs = RVu / ( config->GetSpecificUnits()->GetHbarTimesLambda() * 2 * PI * config->GetSpecificUnits()->GetTimeToHertzFactor()) ; ///Energy in hbar * Hz
	  double IVu = imag(config->GetSpecificUnits()->KValueToEnergy(printVector[id])); ///Energy in custom unit. = -Gamma/2
	  double IVrte = 2.0*abs(IVu) / (config->GetSpecificUnits()->GetHbarTimesLambda() * config->GetSpecificUnits()->GetTimeToHertzFactor()); /// Gamma / hbar 

	  if(imag(printVector[id]) > 1E-5 && abs(arg(printVector[id])-PI/2) < 1E-2 )
		{
		  vPrint(1,"Bound state: k = %+6.10fi [%s]^(-1)    =>   E = %+6.10f %s  (= %6.10f Hz * h) \n", 
				 imag(printVector[id]), 
				 config->GetSpecificUnits()->GetLengthUnitName().c_str(), 
				 RVu, 
				 config->GetSpecificUnits()->GetEnergyUnitName().c_str(), 
				 RVs
				 );
		}
	  else if((imag(printVector[id]) < -1E-6 && arg(printVector[id]) < 0.0 && arg(printVector[id]) > -1.0*PI/2 ))
		{
		  vPrint(1,"Resonant state (index %ld): k = [ %+6.10f %+6.10fi ] [%s]^(-1)    =>    E = [ %+6.10f %+6.10fi ] %s  (= %6.10f Hz * h,  decayRate= %6.10f s^{-1})\n", 
				 idxVector[id],
				 real(printVector[id]), 
				 imag(printVector[id]), 
				 config->GetSpecificUnits()->GetLengthUnitName().c_str(), 
				 RVu, 
				 imag(config->GetSpecificUnits()->KValueToEnergy(printVector[id])), 
				 config->GetSpecificUnits()->GetEnergyUnitName().c_str(), 
				 RVs, 
				 IVrte
				 );
		}
	  else
		{
		  vPrint(1,"Other state (index %ld): k = [ %+6.10f %+6.10fi ] [%s]^(-1)  =>   E = [ %+6.10f + %6.10fi ] %s (= %6.10f Hz * h,  decayRate= %6.10f s^{-1}) \n", 
				 idxVector[id],
				 real(printVector[id]), 
				 imag(printVector[id]),
				 config->GetSpecificUnits()->GetEnergyUnitName().c_str(), 
				 RVu, 
				 IVu, 
				 config->GetSpecificUnits()->GetEnergyUnitName().c_str(), 
				 RVs, 
				 IVrte
				 );
		}
	}

}

void OutputProcessor::WriteInterestingKPointsToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("InterestingPointsFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty InterestingPoints filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving interesting k-points...");

  vector<ComplexDouble> printVector = FindInterestingKPoints();

  FILE * fout = AssuredFopen(fileName);  
  for(ulong i = 0; i<printVector.size(); ++i)
	{
	  fprintf(fout,"%+13.10e %+13.10e\n", real(printVector[i]), imag(printVector[i]));
	}

  fclose(fout);
  fout = NULL;
  vPrint(4, "done\n");
}


vector<ulong> OutputProcessor::FindInterestingKPointIndex(bool forceFilter) const
{
    ///Retrieve k-points corresponding to the basis: those are the "filter", which means that
  ///points too close to them are considered "uninteresting".
  ///Filtering only enabled for 1 particle case however, otherwise they must be explicitly specified in the config file.
  vector<ulong> toReturn;

  if(config->GetNumberOfParticles() == 1 || forceFilter)
	{
	  vector<ComplexDouble> filterVector;
	  for(ulong i = 0; i<config->GetKCurve()->GetNumberOfSegments(); ++i)
		{
		  const vector<pair<ComplexDouble, ComplexDouble> > * rule = config->GetKCurve()->GetSegmentRule(i);
		  for(vector<pair<ComplexDouble, ComplexDouble> >::const_iterator it = rule->begin(); it != rule->end(); ++it)
			{
			  filterVector.push_back(it->first);
			}
		}
	  
	  
	  double imagScale = 1E-8; ///The maximum imaginary part (without sign)
	  for(ulong i = 0; i<filterVector.size(); ++i)
		{
		  imagScale = MAX(imagScale, abs(imag(filterVector[i])));
		}
	  
	  
	  ///Filter out the uninteresting values and return the interesting ones.
	  ulong eigenCounter = 0;
	  for(vector<EigenPair>::const_iterator it = eigenData->EigenPairs.begin(); it!=eigenData->EigenPairs.end(); ++it)
		{
		  ComplexDouble kToPrint = config->GetSpecificUnits()->EnergyToKValue(it->Eigenvalue);
		  
		  ///Now apply filter rule.
		  
		  
		  
		  ///All k-values on the imaginary axis are interesting.
		  if(imag(kToPrint) > 1E-5 && abs(arg(kToPrint)-PI/2) < 1E-2 )
			{
			  toReturn.push_back(eigenCounter);
			}
		  
		  ///Right half plane points might be interesting.
		  else if(real(kToPrint) > 1E-5 )
			{
			  bool tooClose = false;
			  bool tooRight = true; ///Ignore rightmost points: they are never interesting.
			  for(ulong i = 0; i<filterVector.size(); ++i)
				{
				  ComplexDouble d1 = filterVector[i] - kToPrint;
				  ComplexDouble d2 = 0;
				  
				  if(i>0)
					{
					  d2 += filterVector[i]-filterVector[i-1];
					}
				  
				  if(i+1<filterVector.size())
					d2 += filterVector[i+1]-filterVector[i];
				  if(i > 0 && i+1 < filterVector.size())
					d2 /= 2.0;
				  
				  if(i +1 < filterVector.size())
					d2 *= 1.1;
				  
				  
				  if(abs(real(d1)) < abs(real(d2)) && ((abs(imag(d1)) < abs(imag(d2)))  || abs(imag(d1)) < 0.05 * imagScale) )
					{
					  tooClose = true;
					  break;
					}
				  
				  if(real(filterVector[i]) > real(kToPrint))
					{
					  tooRight = false;
					}
				  
				}
			  if(!tooClose && !tooRight)
				toReturn.push_back(eigenCounter);
			}
		  ++eigenCounter;
		}
	  if(eigenCounter != eigenData->EigenPairs.size() )
		throw RLException("Consistency error: invalid number of eigenvalues.");
	}

  vector<ComplexDouble> extraInterestingPoints = config->GetExtraInterestingPoints();

  ///Append some other points that are "enforced interesting".
  for(vector<ComplexDouble>::const_iterator ip = extraInterestingPoints.begin(); ip!=extraInterestingPoints.end(); ++ip)
	{
	  double minDistance = 1E99;
	  int bestMatch = -1;
	  

	  ulong loopCount = 0;
	  for(vector<EigenPair>::const_iterator it = eigenData->EigenPairs.begin(); it!=eigenData->EigenPairs.end(); ++it)
		{
		  ComplexDouble kToPrint = config->GetSpecificUnits()->EnergyToKValue(it->Eigenvalue);
		  double locDist = abs(kToPrint-*ip);
		  if(minDistance > locDist)
			{
			  minDistance = locDist;
			  bestMatch = loopCount;
			}
		  ++loopCount;
		}
	  if(DBL_EQUAL(bestMatch, -1))
		throw RLException("Best match was never found. This should never happen.");
	  toReturn.push_back((ulong)bestMatch);
	}
  return toReturn;
}




vector<ComplexDouble> OutputProcessor::FindInterestingKPoints(bool forceFilter) const
{
  vector<ComplexDouble> toReturn;
  vector<ulong> ikpIndex = FindInterestingKPointIndex(forceFilter);
  for(vector<ulong>::const_iterator it = ikpIndex.begin(); it!=ikpIndex.end(); ++it)
	{
	  if(*it >= eigenData->EigenPairs.size())
		throw RLException("Consistency error: interesting index was larger than eigenvalue vector size.");
	  toReturn.push_back(config->GetSpecificUnits()->EnergyToKValue(eigenData->EigenPairs[*it].Eigenvalue));
	}

  return toReturn;
}





void OutputProcessor::WriteInterestingOneParticleWavefunctionsToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("WavefunctionsFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty Wavefunctions filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving wavefunctions...\n");

  ///Some setup.
  vector<ulong> interestingIndexes = FindInterestingKPointIndex();
  vector<ComplexDouble> interestingVector = FindInterestingKPoints();


  CompositeBasisFunction * myCompositeBasisFunction = NULL;


  if(! config->GetHarmonicOverride())
	{
	  myCompositeBasisFunction = new CompositeBasisFunction(config->GetBasisFunctions(), 
															eigenData, 
															config->GetKCurve());
	}
  else
	{
	  myCompositeBasisFunction = new CompositeBasisFunction(config->GetHarmonicBasisFunction(), 
															eigenData,
															0
															);
	}
  
  vector<double> xValues;
  for(double x = config->GetMinWavefunctionX(); x <= config->GetMaxWavefunctionX(); x += config->GetWavefunctionStepsizeX())
	xValues.push_back(x);

  ///Compute the wavefunctions : this is where the magic is.
  vector<vector<ComplexDouble> > wavefunctionValues(interestingIndexes.size(), vector<ComplexDouble>(xValues.size(), 0.0));
  for(ulong iPtr = 0; iPtr < interestingVector.size(); ++iPtr)
	{
	  for(ulong i = 0; i<xValues.size(); ++i)
		{
		  wavefunctionValues.at(iPtr).at(i) += myCompositeBasisFunction->Eval(xValues.at(i), interestingIndexes.at(iPtr));
		}
	  
	  //As a verification, compute the integral of the wavefunction and print it.
	  double sqSum = SquareSum(wavefunctionValues.at(iPtr));
	  vPrint(2, "k = %+2.5f%+2.5f => Graph area: %+13.10e\n", real(interestingVector.at(iPtr)), imag(interestingVector.at(iPtr)), sqSum * config->GetWavefunctionStepsizeX());
	}

  ///Save to file.
  FILE * fout = AssuredFopen(fileName);

  fprintf(fout, "#x-value");
  for(ulong j = 0; j<wavefunctionValues.size(); ++j)
	{
	  fprintf(fout, " Real_wave_%ld(k=%+2.5f%+2.5fi) Imag_wave_%ld", j, real(interestingVector[j]), imag(interestingVector[j]),j);
	}
  fprintf(fout, "\n");
  for(ulong i = 0; i<xValues.size(); ++i)
	{
	  fprintf(fout, "%+13.10e", xValues.at(i));
	  for(ulong j = 0; j < wavefunctionValues.size(); ++j)
		{
		  fprintf(fout, " %+13.10e %+13.10e", 
				  real(wavefunctionValues.at(j).at(i)), 
				  imag(wavefunctionValues.at(j).at(i))
				  );
		}
	  fprintf(fout, "\n");
	}
  fclose(fout); fout = NULL;

  vPrint(4, "Saving wavefunctions...done\n");
}


double OutputProcessor::SquareSum( const vector<ComplexDouble> & toSum)
{
  double sqSum = 0;
  for(vector<ComplexDouble>::const_iterator ip = toSum.begin(); ip!=toSum.end(); ++ip)
	{
	  sqSum += pow(abs(*ip), 2.0);
	}
  return sqSum;
}

bool OutputProcessor::SaveMatrix(CMatrix * toSave) const
{
  string fileName = config->GetOutputFilenames()->Get("MatrixFile");


  if(fileName.empty() )
	{
	  vPrint(4, "Empty Matrix filename, not saving.\n");
	  return false;
	}


  vPrint(4, "Saving matrix...");
  FILE * fout = AssuredFopen(fileName);

  for(ulong i = 0; i<toSave->Rows(); ++i)
	{
	  for(ulong j = 0; j<toSave->Columns(); ++j)
		{
		  fprintf(fout, "%ld %ld %+13.10e %+13.10e\n", i, j, real(toSave->Element(i, j)), imag(toSave->Element(i, j)));
		}
	}
  fclose(fout);

  vPrint(4, "done\n");
  return true;
}



void OutputProcessor::WriteInterestingRelativeTwoParticleWavefunctionsToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("WavefunctionsFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty Wavefunctions filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving relative wavefunctions...\n");

  if(myCompositeBasisFunctions == NULL)
	{
	  throw RLException("Must set CompositeBasisFunctions for OutputProcessor before writing to file.");
	}
  if(myCompositeBasisFunctions->size() < 2)
	{
	  vPrint(2, "Too few composite basis functions, not saving.\n");
	  return;
	}

  ///Some setup.
  vector<ulong> interestingIndexes = FindInterestingKPointIndex();
  vector<ComplexDouble> interestingVector = FindInterestingKPoints();


  vector<double> xrelValues;

  for(double x = config->GetMinWavefunctionX(); x <= config->GetMaxWavefunctionX(); x += config->GetWavefunctionStepsizeX())
    xrelValues.push_back(x);


  vector<pair<double, double> > wflRule = 
	LegendreRule::GetRule(xrelValues.size(), 
						  config->GetMinWavefunctionX(),
	                      config->GetMaxWavefunctionX());



  ///Compute the wavefunctions : this is where the magic is.
  vector<vector<double> > wavefunctionValues(interestingIndexes.size(), vector<double>(xrelValues.size(), 0.0));

  ulong N1 = myCompositeBasisFunctions->at(0)->GetSize();
  ulong N2 = myCompositeBasisFunctions->at(1)->GetSize();


  for(ulong i = 0; i<interestingIndexes.size(); ++i)
  {
	vPrint(4, "Wavefunction # %ld (of %ld)\n", i, interestingIndexes.size());
	vPrint(4, "Progress: 00.00%%");

	for(ulong j = 0; j<xrelValues.size(); ++j)
		{
		  vPrint(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
		  vPrint(4, "Progress: %04.02f%%", (double)100*((double)j/xrelValues.size()));
		  double wval = 0;
		  for(ulong xabsP = 0; xabsP < wflRule.size(); ++xabsP)
			{
			  ComplexDouble sum = 0.0;			  
			  double xabs = wflRule.at(xabsP).first;
			  double weight = wflRule.at(xabsP).second;

			  for(ulong k = 0; k<eigenData->EigenPairs.at(interestingIndexes.at(i)).Eigenvector.size(); ++k)
				{
				  double xrel = xrelValues.at(j);

				  
				  ulong a = k / N1;
				  ulong b = k % N2;
				  
				  sum += 
					eigenData->EigenPairs.at(interestingIndexes.at(i)).Eigenvector.at(k) *
					myCompositeBasisFunctions->at(0)->Eval((xabs-xrel)/sqrt(2.0), a) *
					myCompositeBasisFunctions->at(1)->Eval((xabs+xrel)/sqrt(2.0), b);
				}
			  wval += pow(abs(sum), 2.0) * weight;
			}
		  wavefunctionValues.at(i).at(j) += wval;
		}
	double totalSum = 0;
	for(ulong a = 0; a<xrelValues.size(); ++a)
	  {
		totalSum += wavefunctionValues.at(i).at(a);
	  }
	vPrint(2, "Total rel-graph area: %+13.6e\n", totalSum * (config->GetMaxWavefunctionX() - config->GetMinWavefunctionX()) / xrelValues.size());
	}


  ///Save to file.
  FILE * fout = AssuredFopen(fileName);
  
  fprintf(fout, "#x-value");
  for(ulong j = 0; j<wavefunctionValues.size(); ++j)
	{
	  fprintf(fout, " WaveABS_%ld(k=%+2.5f%+2.5fi) zero-padding", j, real(interestingVector[j]), imag(interestingVector[j]));
	}
  fprintf(fout, "\n");
  for(ulong i = 0; i<xrelValues.size(); ++i)
	{
	  fprintf(fout, "%+13.10e", xrelValues[i]);
	  for(ulong j = 0; j < wavefunctionValues.size(); ++j)
		{
		  fprintf(fout, " %+13.10e 0.0000000", sqrt(wavefunctionValues.at(j).at(i)));
		}
	  fprintf(fout, "\n");
	}
  fclose(fout);
  fout = NULL;


  vPrint(4, "Saving relative wavefunctions...done\n");
}




void OutputProcessor::WriteProductTwoParticleWavefunctionToFile() const
{

  string fileName = config->GetOutputFilenames()->Get("ProductWavefunction");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty ProductWavefunction filename, not saving.\n");
	  return;
	}

  if(myCompositeBasisFunctions->size() < 2)
	{
	  vPrint(2, "Too few composite basis functions, not saving.\n");
	  return;
	}



  FILE * fout = AssuredFopen(fileName);

  vPrint(4, "Saving most interesting product wavefunction...\n");

  vector<double> xrelValues;
  
  for(double x = config->GetMinWavefunctionX(); x <= config->GetMaxWavefunctionX(); x += config->GetWavefunctionStepsizeX())
    xrelValues.push_back(x);

  vector<vector<double> > mesh(xrelValues.size(), vector<double>(xrelValues.size(), 0.0));
  
  vector<ulong> interestingIndexes = FindInterestingKPointIndex();
  if(interestingIndexes.empty())
	return;
 
  
  ulong N1 = myCompositeBasisFunctions->at(0)->GetSize();
  ulong N2 = myCompositeBasisFunctions->at(1)->GetSize();

  if(eigenData->EigenPairs.size() != N1 * N2)
	{
	  throw RLException("Size mismatch: N1 = %ld, N2 = %ld but there are %ld eigenvectors.", N1, N2, eigenData->EigenPairs.size());
	}
  
  double area = 0;
  vPrint(4, "Progress: 00.00%%");
  for(ulong ux1 = 0; ux1 < xrelValues.size(); ++ux1)
	{
	  vPrint(4, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\033[K");
	  vPrint(4, "Progress: %04.02f%%", (double)100*((double)ux1/xrelValues.size()));
	  for(ulong ux2 = 0; ux2 < xrelValues.size(); ++ux2)
		{
		  for(ulong i = 0; i<eigenData->EigenPairs.size(); ++i)
			{
			  ulong a = i / N1;
			  ulong b = i % N2;
			  double x1 = xrelValues[ux1];
			  double x2 = xrelValues[ux2];
			  mesh.at(ux1).at(ux2) += real(
										     myCompositeBasisFunctions->at(0)->Eval(x1, a) 
										   * myCompositeBasisFunctions->at(1)->Eval(x2, b) 
										   * eigenData->EigenPairs.at(interestingIndexes.front()).Eigenvector.at(i)
										   );
			  
			}
		  area += pow(abs(mesh.at(ux1).at(ux2)), 2.0) * pow(abs(config->GetWavefunctionStepsizeX()), 2.0);
		}
	}
  vPrint(2, "2D probability area: %+13.10e\n", area);

  
  for(ulong a = 0; a<xrelValues.size(); ++a)
	{
	  for(ulong b = 0; b<xrelValues.size(); ++b)
		{
		  fprintf(fout, "%13.10e ", mesh[a][b]);
		}
	  fprintf(fout, "\n");
	}
  
  
  fclose(fout); fout = NULL;

  vPrint(4, "Saving most interesting product wavefunction...done\n");

}
