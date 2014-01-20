#include "ExternalLauncher.hh"

ExternalLauncher::ExternalLauncher(ComputeConfig * _config)
  :config(_config)
{
  
}

void ExternalLauncher::Launch() const
{
  if(config->GetAutoPlotPotential() && config->GetOutputFilenames()->Get("PotentialFile").length() > 1)
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
  if(config->GetAutoPlotKCurve() && config->GetOutputFilenames()->Get("KCurveFile").length() > 1)
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./KPlot.sh \"%s\" \"%s\" \"%s\"", 
			  config->GetOutputFilenames()->Get("KCurveFile").c_str(),
			  config->GetOutputFilenames()->Get("KFoundFile").c_str(),
			  config->GetSpecificUnits()->GetLengthUnitName().c_str()
			  );

	  system(buffer);
	}

  if(config->GetAutoPlotWavefunctions() && config->GetOutputFilenames()->Get("WavefunctionsFile").length() > 1)
	{
	  char buffer[4000];
	  sprintf(buffer, "gnuplot/./WFplot.sh \"%s\" \"%s\" ", 
			  config->GetOutputFilenames()->Get("WavefunctionsFile").c_str(),
			  config->GetSpecificUnits()->GetLengthUnitName().c_str()
			  );
	  
	  system(buffer);
	}
}
