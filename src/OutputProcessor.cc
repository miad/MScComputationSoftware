#include "OutputProcessor.hh"

OutputProcessor::OutputProcessor(ComputeConfig * _config)
  :config(_config), eigenData(NULL)
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


void OutputProcessor::WritePostOutput() const
{
  if(eigenData == NULL)
	throw RLException("OutputProcessor error: eigenData was null.");

  WritePotentialToFile();
  WritePotentialPrecisionToFile();
  WriteKCurveToFile();
  WriteKFoundToFile();
  WriteInterestingKPointsVerbosely();
  if(config->GetNumberOfParticles() == 1)
	{
	  WriteInterestingKPointsToFile();
	  WriteInterestingOneParticleWavefunctionsToFile();
	}
  else if(config->GetNumberOfParticles() == 2)
	{
	  WriteInterestingTwoParticleWavefunctionsToFile();
	}
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

  for(uint i = 0; i<toPrint->GetNumberOfSegments(); ++i)
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

  if(fileName.empty() )
	{
	  vPrint(4, "Empty KFound filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving k-values for found eigenvalues...");

  FILE * fout = AssuredFopen(fileName);

  for(vector<ComplexDouble>::const_iterator it = eigenData->Eigenvalues.begin(); it!=eigenData->Eigenvalues.end(); ++it)
	{
	  ComplexDouble kToPrint = EnergyToKValue(*it);

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

ComplexDouble OutputProcessor::EnergyToKValue(const ComplexDouble & energy) const
{
  ///If numerical stability is mean to us, then rotate.
  ComplexDouble toReturn = sqrt(energy * 2. * (config->GetSpecificUnits()->GetMassOverLambda2())) / 
	(config->GetSpecificUnits()->GetHbarTimesLambda());
  
  if( (abs(imag(toReturn)) > 1E2*abs(real(toReturn))  && imag(toReturn) < 0))
	{
	  toReturn *= -1.0;
	}

  return toReturn;
}

ComplexDouble OutputProcessor::KValueToEnergy(const ComplexDouble & kValue) const
{
  return pow(kValue * (config->GetSpecificUnits()->GetHbarTimesLambda()), 2.0) /
	(2.0 * (config->GetSpecificUnits()->GetMassOverLambda2()) );
}




void OutputProcessor::WriteInterestingKPointsVerbosely() const
{
  vector<ComplexDouble> printVector = FindInterestingKPoints();

  for(vector<ComplexDouble>::const_iterator it = printVector.begin(); it!=printVector.end(); ++it)
	{
	  if(imag(*it) > 1E-5 && abs(arg(*it)-PI/2) < 1E-2 )
		{
		  vPrint(1,"Bound state: k = %+6.10lfi [%s]^(-1)    =>   E = %+6.10lf %s \n", imag(*it), config->GetSpecificUnits()->GetLengthUnitName().c_str(), real(KValueToEnergy(*it)), config->GetSpecificUnits()->GetEnergyUnitName().c_str());
		}
	  else if((imag(*it) < -1E-6 && arg(*it) < 0.0 && arg(*it) > -1.0*PI/2 ))
		{
		  vPrint(1,"Resonant state: k = [ %+6.10lf %+6.10lfi ] [%s]^(-1)    =>    E = [ %+6.10lf %+6.10lfi] %s\n", real(*it), imag(*it), config->GetSpecificUnits()->GetLengthUnitName().c_str(), real(KValueToEnergy(*it)), imag(KValueToEnergy(*it)), config->GetSpecificUnits()->GetEnergyUnitName().c_str());
		}
	  else
		{
		  vPrint(1,"Other state: k = %+6.10lf %+6.10lfi\n", real(*it), imag(*it));
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
  vector<uint> interestingIndexes = FindInterestingKPointIndex();

  FILE * fout = AssuredFopen(fileName);
  
  
  if(printVector.size() != interestingIndexes.size())
	{
	  throw RLException("There were %d print vectors but %d interestingIndexes, they should match.\n", printVector.size(), interestingIndexes.size());
	}
  
  for(uint i = 0; i<printVector.size(); ++i)
	{
	  fprintf(fout,"%+13.10e %+13.10e", real(printVector[i]), imag(printVector[i]));
	  double sum = -1;
	  vector<double> ratio = GetOneParticleBasisRatio(interestingIndexes[i], sum);
	  for(vector<double>::const_iterator it = ratio.begin(); it!=ratio.end(); ++it)
		{
		  fprintf(fout, " %+13.10e", *it);
		}
	  fprintf(fout, " %+13.10e\n", sum);
	}

  
  fclose(fout);
  fout = NULL;
  vPrint(4, "done\n");
}


vector<uint> OutputProcessor::FindInterestingKPointIndex() const
{
    ///Retrieve k-points corresponding to the basis: those are the "filter", which means that
  ///points too close to them are considered "uninteresting".
  vector<ComplexDouble> filterVector;
  for(uint i = 0; i<config->GetKCurve()->GetNumberOfSegments(); ++i)
	{
	  const vector<pair<ComplexDouble, ComplexDouble> > * rule = config->GetKCurve()->GetSegmentRule(i);
	  for(vector<pair<ComplexDouble, ComplexDouble> >::const_iterator it = rule->begin(); it != rule->end(); ++it)
		{
		  filterVector.push_back(it->first);
		}
	}
  
  ///Filter out the uninteresting values and return the interesting ones.
  vector<uint> toReturn;
  uint eigenCounter = 0;
  for(vector<ComplexDouble>::const_iterator it = eigenData->Eigenvalues.begin(); it!=eigenData->Eigenvalues.end(); ++it)
	{
	  ComplexDouble kToPrint = EnergyToKValue(*it);

	  ///Now apply filter rule.
	  if(imag(kToPrint) > 1E-5 && abs(arg(kToPrint)-PI/2) < 1E-2 )
		toReturn.push_back(eigenCounter);
	  else if(real(kToPrint) > 1E-5 )
		{
		  bool tooClose = false;
		  for(uint i = 0; i<filterVector.size(); ++i)
			{
			  double d1 = abs(filterVector[i] - kToPrint);
			  double d2 = 0;
			  if(i>0)
				d2 += abs(filterVector[i-1]-filterVector[i]);
			  if(i+1<filterVector.size())
				d2 += abs(filterVector[i+1]-filterVector[i]);
			  if(i > 0 && i+1 < filterVector.size())
				d2 /= 2;
			  if(i +1 < filterVector.size())
				d2 *= 0.6;

			  if(d1 < d2)
				{
				  tooClose = true;
				  break;
				}
			}
		  if(!tooClose)
			toReturn.push_back(eigenCounter);
		}
	  ++eigenCounter;
	}
  if(eigenCounter != eigenData->Eigenvalues.size() )
	throw RLException("Consistency error: invalid number of eigenvalues.");
  
  vector<ComplexDouble> extraInterestingPoints = config->GetExtraInterestingPoints();

  ///Append some other points that are "enforced interesting".
  for(vector<ComplexDouble>::const_iterator ip = extraInterestingPoints.begin(); ip!=extraInterestingPoints.end(); ++ip)
	{
	  double minDistance = 1E99;
	  uint bestMatch = -1;
	  

	  uint loopCount = 0;
	  for(vector<ComplexDouble>::const_iterator it = eigenData->Eigenvalues.begin(); it!=eigenData->Eigenvalues.end(); ++it)
		{
		  ComplexDouble kToPrint = EnergyToKValue(*it);
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
	  toReturn.push_back(bestMatch);
	}
  return toReturn;
}




vector<ComplexDouble> OutputProcessor::FindInterestingKPoints() const
{
  vector<ComplexDouble> toReturn;
  vector<uint> ikpIndex = FindInterestingKPointIndex();
  for(vector<uint>::const_iterator it = ikpIndex.begin(); it!=ikpIndex.end(); ++it)
	{
	  if(*it >= eigenData->Eigenvalues.size())
		throw RLException("Consistency error: interesting index was larger than eigenvalue vector size.");
	  toReturn.push_back(EnergyToKValue(eigenData->Eigenvalues[*it]));
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
  vector<uint> interestingIndexes = FindInterestingKPointIndex();
  vector<ComplexDouble> interestingVector = FindInterestingKPoints();

  vector<BasisFunction> myBasisFunctions = config->GetBasisFunctions();
  uint numberOfGLPoints = eigenData->Eigenvalues.size() / config->GetBasisFunctions().size();
  if(eigenData->Eigenvalues.size() % config->GetBasisFunctions().size() != 0)
	{
	  throw RLException("Eigenvalue number was not divisible by basis number size.");
	}


  ///Compute the wavefunctions : this is where the magic is.
  vector<vector<ComplexDouble> > wavefunctionValues;
  for(vector<uint>::const_iterator it = interestingIndexes.begin(); it!=interestingIndexes.end(); ++it)
	{
	  wavefunctionValues.push_back(vector<ComplexDouble>());
	  vector<ComplexDouble> eigVect = GetReshapedEigenvector(*it);
	  for(double x = config->GetMinWavefunctionX(); x <= config->GetMaxWavefunctionX(); x += config->GetWavefunctionStepsizeX())
		{
		  wavefunctionValues.back().push_back(ComplexDouble(0.0,0.0));
		  for(uint i = 0; i<eigVect.size(); ++i)
			{
			  uint curvePointer = i % numberOfGLPoints;
			  uint basisPointer = i / numberOfGLPoints;
			  uint curveSegment = config->GetKCurve()->SegmentIndexFromGLNumber(curvePointer);
			  ComplexDouble kVal = config->GetKCurve()->GetRuleValue(curveSegment, curvePointer);
			  ComplexDouble wK = config->GetKCurve()->GetRuleWeight(curveSegment, curvePointer);
			  
			  wavefunctionValues.back().back() += 
				(eigVect[i] * sqrt(wK) * 
				 myBasisFunctions[basisPointer].Eval(x, kVal)
				 );
			}
		}

	  //As a verification, compute the integral of the wavefunction and print it.
	  double sqSum = 0;
	  for(vector<ComplexDouble>::const_iterator ip = wavefunctionValues.back().begin(); ip!=wavefunctionValues.back().end(); ++ip)
		{
		  sqSum += pow(abs(*ip), 2) * config->GetWavefunctionStepsizeX();
		}
	  vPrint(2, "k = %+2.5lf%+2.5lf => Graph area: %+13.10e\n", real(EnergyToKValue(eigenData->Eigenvalues[*it])), imag(EnergyToKValue(eigenData->Eigenvalues[*it])), sqSum);

	}


  ///Save to file.
  FILE * fout = AssuredFopen(fileName);
  
  fprintf(fout, "#x-value");
  for(uint j = 0; j<wavefunctionValues.size(); ++j)
	{
	  fprintf(fout, " Real_wave_%d(k=%+2.5lf%+2.5lfi) Imag_wave_%d", j, real(interestingVector[j]), imag(interestingVector[j]),j);
	}
  fprintf(fout, "\n");
  uint i = -1;
  for(double x = config->GetMinWavefunctionX(); x <= config->GetMaxWavefunctionX(); x += config->GetWavefunctionStepsizeX())
	{
	  ++i;
	  fprintf(fout, "%+13.10e", x);
	  for(uint j = 0; j < wavefunctionValues.size(); ++j)
		{
		  fprintf(fout, " %+13.10e %+13.10e", real(wavefunctionValues.at(j).at(i)), imag(wavefunctionValues.at(j).at(i)));
		}
	  fprintf(fout, "\n");
	}
  fclose(fout);
  fout = NULL;


  ///Compute basis element decomposition and print it.
  int loopCounter = 0;
  for(vector<uint>::const_iterator it = interestingIndexes.begin(); it!=interestingIndexes.end(); ++it)
	{
	  double totalSum = -1;
	  vector<double> basisRatio = GetOneParticleBasisRatio(*it, totalSum);

	  vPrint(3, "Basis element decomposition for k = %+10.6lf%+10.6lfi:\n", real(interestingVector[loopCounter]), imag(interestingVector[loopCounter]));
	  ++loopCounter;
	  
	  for(uint i = 0; i<basisRatio.size(); ++i)
		{
		  vPrint(3, "\t%s : \t%3.2lf%%\n", myBasisFunctions[i].GetName(), basisRatio[i]*100);
		}
	  vPrint(3, "Element square sum: %+13.10e\n\n", totalSum);
	}
  
  vPrint(4, "Saving wavefunctions...done\n");
}


vector<double> OutputProcessor::GetOneParticleBasisRatio(uint eigenIndex, double & totalSum) const
{
  totalSum = 0;
  if(eigenIndex > eigenData->Eigenvectors.size())
	{
	  throw RLException("Too large eigenIndex used.");
	}

  uint numberOfBasisVectors = config->GetBasisFunctions().size();

  uint numberOfGLPoints = config->GetKCurve()->GetTotalNumberOfGLPoints();


  if(eigenData->Eigenvectors[eigenIndex].size() % numberOfBasisVectors != 0)
	{
	  throw RLException("Vector length and # basis vectors does not match.");
	}

  vector<double> basisRatio;
  for(uint i = 0; i<numberOfBasisVectors; ++i)
	{
	  basisRatio.push_back(0);
	}
  for(uint i = 0; i<eigenData->Eigenvectors[eigenIndex].size(); ++i)
	{
	  uint curvePointer = i % numberOfGLPoints;
	  uint basisPointer = i / numberOfGLPoints;
	  uint curveSegment = config->GetKCurve()->SegmentIndexFromGLNumber(curvePointer);
	  basisRatio.at(basisPointer) += 
		abs(config->GetKCurve()->GetRuleWeight(curveSegment, curvePointer))*
		pow(abs(eigenData->Eigenvectors.at(eigenIndex).at(i) * config->GetKCurve()->GetRuleWeight(curveSegment, curvePointer)),2);
	}
  

  
  for(uint i = 0; i<basisRatio.size(); ++i)
	{
	  totalSum += basisRatio.at(i);
	}
  for(uint i = 0; i<basisRatio.size(); ++i)
	{
	  basisRatio.at(i) /= totalSum;
	}
  return basisRatio;
}




vector<ComplexDouble> OutputProcessor::GetReshapedEigenvector(uint index) const
{
  if(index >= eigenData->Eigenvectors.size())
	{
	  throw RLException("Tried to access eigenvector %d, but max is %d.", index, eigenData->Eigenvectors.size());
	}
  vector<ComplexDouble> toReturn = eigenData->Eigenvectors[index];

  ComplexDouble normalizationFactor = ComplexDouble(0.0,0.0);
  for(uint i = 0; i<toReturn.size(); ++i)
	{
 	  normalizationFactor += pow(toReturn[i], 2.0);
	}
  ComplexDouble toDivideBy = sqrt(normalizationFactor);
  for(vector<ComplexDouble>::iterator it = toReturn.begin(); it!=toReturn.end(); ++it)
	{
	  *it /= toDivideBy;
	}
  return toReturn;
}


void OutputProcessor::SaveMatrix(CMatrix * toSave) const
{
  string fileName = config->GetOutputFilenames()->Get("MatrixFile");


  if(fileName.empty() )
	{
	  vPrint(4, "Empty Matrix filename, not saving.\n");
	  return;
	}


  vPrint(4, "Saving matrix...");
  FILE * fout = AssuredFopen(fileName);

  for(uint i = 0; i<toSave->Rows(); ++i)
	{
	  for(uint j = 0; j<toSave->Columns(); ++j)
		{
		  fprintf(fout, "%d %d %+13.10e %+13.10e\n", i, j, real(toSave->Element(i, j)), imag(toSave->Element(i, j)));
		}
	}
  fclose(fout);

  vPrint(4, "done\n");
}





void OutputProcessor::WriteInterestingTwoParticleWavefunctionsToFile() const
{
  string fileName = config->GetOutputFilenames()->Get("WavefunctionsFile");

  if(fileName.empty() )
	{
	  vPrint(4, "Empty Wavefunctions filename, not saving.\n");
	  return;
	}

  vPrint(4, "Saving wavefunctions...\n");

  ///Some setup.
  vector<uint> interestingIndexes = FindInterestingKPointIndex();
  vector<ComplexDouble> interestingVector = FindInterestingKPoints();

  vector<BasisFunction> myBasisFunctions = config->GetBasisFunctions();
  uint numberOfGLPoints = config->GetKCurve()->GetTotalNumberOfGLPoints();

  uint N = myBasisFunctions.size() * numberOfGLPoints; //Subspace dimension.


  if(eigenData->Eigenvalues.size() % N != 0)
	{
	  throw RLException("Eigenvalue number was invalid.");
	}


  ///Compute the wavefunctions : this is where the magic is.
  vector<vector<vector<ComplexDouble> > > wavefunctionValues;
  for(vector<uint>::const_iterator it = interestingIndexes.begin(); it!=interestingIndexes.end(); ++it)
	{
	  wavefunctionValues.push_back(vector<vector<ComplexDouble> >());

	  vector<ComplexDouble> eigVect = GetReshapedEigenvector(*it);

	  for(double x1 = config->GetMinWavefunctionX(); x1 <= config->GetMaxWavefunctionX(); x1 += config->GetWavefunctionStepsizeX())
		{
		  wavefunctionValues.back().push_back(vector<ComplexDouble>() );
		  for(double x2 = config->GetMinWavefunctionX(); x2 <= config->GetMaxWavefunctionX(); x2 += config->GetWavefunctionStepsizeX())
			{
			  wavefunctionValues.back().back().push_back(ComplexDouble(0.0,0.0));
			  for(uint i = 0; i<eigVect.size(); ++i)
				{


				  uint a = i/N;
				  uint b = i % N;
				  
				  uint curvePointerA = a % numberOfGLPoints;
				  uint basisPointerA = a / numberOfGLPoints;
				  uint curveSegmentA = config->GetKCurve()->SegmentIndexFromGLNumber(curvePointerA);
				  
				  ComplexDouble kA = config->GetKCurve()->GetRuleValue(curveSegmentA, curvePointerA);
				  ComplexDouble wA = config->GetKCurve()->GetRuleWeight(curveSegmentA, curvePointerA);

				  uint curvePointerB = b % numberOfGLPoints;
				  uint basisPointerB = b / numberOfGLPoints;
				  uint curveSegmentB = config->GetKCurve()->SegmentIndexFromGLNumber(curvePointerB);
				  
				  ComplexDouble kB = config->GetKCurve()->GetRuleValue(curveSegmentB, curvePointerB);
				  ComplexDouble wB = config->GetKCurve()->GetRuleWeight(curveSegmentB, curvePointerB);
				  
				  wavefunctionValues.back().back().back() += 
					eigVect.at(i) * sqrt(wA)* sqrt(wB) * 
					myBasisFunctions.at(basisPointerA).Eval(x1, kA) * myBasisFunctions.at(basisPointerB).Eval(x2, kB);

				}
			}
		}
	  




	  //As a verification, compute the integral of the wavefunction and print it.
	  double sqSum = 0;
	  for(vector<vector<ComplexDouble> >::const_iterator ip = wavefunctionValues.back().begin(); ip!=wavefunctionValues.back().end(); ++ip)
		{
		  for(vector<ComplexDouble>::const_iterator iq = ip->begin(); iq != ip->end(); ++iq)
			{
			  sqSum += pow(abs(*iq), 2) * pow(config->GetWavefunctionStepsizeX(),2);
			}
		}
	  vPrint(2, "k = %+2.5lf%+2.5lf => Graph area: %+13.10e\n", real(EnergyToKValue(eigenData->Eigenvalues[*it])), imag(EnergyToKValue(eigenData->Eigenvalues[*it])), sqSum);

	}


  ///Save to file.
  FILE * fout = AssuredFopen(fileName);
  
  fprintf(fout, "#x1 #x2");
  for(uint j = 0; j<wavefunctionValues.size(); ++j)
	{
	  fprintf(fout, " Abs2_wave_%d(k=%+2.5lf%+2.5lfi)", j, real(interestingVector[j]), imag(interestingVector[j]));
	}
  fprintf(fout, "\n");
  uint i = -1;
  for(double x1 = config->GetMinWavefunctionX(); x1 <= config->GetMaxWavefunctionX(); x1 += config->GetWavefunctionStepsizeX())
	{
	  ++i;
	  uint j = -1;
	  for(double x2 = config->GetMinWavefunctionX(); x2 <= config->GetMaxWavefunctionX(); x2 += config->GetWavefunctionStepsizeX())
		{
		  ++j;
		  fprintf(fout, "%+13.10e %13.10e", x1, x2);
		  for(uint k = 0; k < wavefunctionValues.size(); ++k)
			{
			  //			  printf("%d %d %d\n", k, i, j);
			  //			  printf("sz %ld\n", wavefunctionValues.at(k).at(i).size());
 			  fprintf(fout, " %+13.10e", pow(abs(wavefunctionValues.at(k).at(i).at(j)),2));
			}
		  fprintf(fout, "\n");
		}
	}
  fclose(fout);
  fout = NULL;
  vPrint(4, "Saving wavefunctions...done\n");
}


