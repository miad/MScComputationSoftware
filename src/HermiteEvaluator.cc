#include "HermiteEvaluator.hh"

bool HermiteEvaluator::autoResize = false;

vector<vector<long double> > HermiteEvaluator::polyCoefficients;

void HermiteEvaluator::Init(uint n, bool _autoResize)
{
  autoResize = _autoResize;
  FillTable(n);
}

long double HermiteEvaluator::HermiteH(uint n, long double x)
{
  if(n >= polyCoefficients.size())
	{
	  if(autoResize)
		{
		  FillTable(n);
		}
	  else
		{
		  throw RLException("No precomputed Hermite polynomial value available for n=%d.", n);
		}
	}
  long double toReturn = 0;
  for(vector<long double>::const_reverse_iterator it = polyCoefficients.at(n).rbegin(); it!= polyCoefficients.at(n).rend(); ++it)
	{
	  toReturn *= x;
	  toReturn += *it;
	}
  return toReturn;
}

void HermiteEvaluator::DeInit()
{
  polyCoefficients.clear();
}

void HermiteEvaluator::FillTable(uint nMax)
{
  /*
	Recursion relation for coefficients:
	a(n+1, k) = 2 a(n, k-1) - 2 n a(n-1, k)          (k>0)
	a(n+1, k) =             - 2 n a(n-1, k)          (k=0)
   */

  for(uint n = polyCoefficients.size(); n<=nMax; ++n)
	{
	  polyCoefficients.push_back(vector<long double>());
	  if(n==0)
		{
		  polyCoefficients[n].push_back(1);
		  continue;
		}
	  else if(n==1)
		{
		  polyCoefficients[n].push_back(0);
		  polyCoefficients[n].push_back(2);
		  continue;
		}
	  for(uint k = 0; k<=n; ++k)
		{
		  long double coeff = 0;
		  coeff += (k>0)?(2.0*polyCoefficients[n-1][k-1]):0;
		  coeff -= (k+1<n)?((n-1)*2.0*polyCoefficients[n-2][k]):0;
		  polyCoefficients[n].push_back(coeff);
		}
	}
}

void HermiteEvaluator::DumpTable()
{
  cout << endl;
  for(uint i = 0; i<polyCoefficients.size(); ++i)
	{
	  for(uint j = 0; j<polyCoefficients[i].size(); ++j)
		{
		  cout << polyCoefficients[i][j] << " ";
		}
	  cout << endl;
	}
}
