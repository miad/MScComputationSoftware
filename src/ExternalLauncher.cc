#include "ExternalLauncher.hh"

ExternalLauncher::ExternalLauncher(ComputeConfig * _config)
  :config(_config)
{
  
}

void ExternalLauncher::Launch() const
{
  if(config->GetAutoPlotPotential())
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./PotentialPlot.sh \"%s\" \"%s\" \"%s\" \"%s\"", 
			  config->GetOutputFilenames()->Get("PotentialFile").c_str(),
			  config->GetOutputFilenames()->Get("PotentialPrecisionFile").c_str(),
			  config->GetSpecificUnits()->GetLengthUnitName().c_str(),
			  config->GetSpecificUnits()->GetEnergyUnitName().c_str()
			  );
	  system(buffer);
	}
  if(config->GetAutoPlotKCurve())
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./KPlot.sh \"%s\" \"%s\" \"%s\"", 
			  config->GetOutputFilenames()->Get("KCurveFile").c_str(),
			  config->GetOutputFilenames()->Get("KFoundFile").c_str(),
			  config->GetSpecificUnits()->GetLengthUnitName().c_str()
			  );

	  system(buffer);
	}

  if(config->GetAutoPlotWavefunctions())
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./WFplot.sh \"%s\" \"%s\" ", 
			  config->GetOutputFilenames()->Get("WavefunctionsFile").c_str(),
			  config->GetSpecificUnits()->GetLengthUnitName().c_str()
			  );
	  
	  system(buffer);
	}
}
