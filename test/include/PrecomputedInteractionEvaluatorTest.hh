#ifndef PrecomputedInteractionEvaluatorTest_hh
#define PrecomputedInteractionEvaluatorTest_hh 1

#include <iostream>
#include <assert.h>
#include <map>
#include <vector>
#include <cmath>
#include <string.h>
#include <stdio.h>

#include "Globals.hpp"
#include "GenericUnitTest.hh"
#include "Matrix.hpp"
#include "PrecomputedInteractionEvaluator.hh"
#

#include <algorithm>

using namespace std;


class PrecomputedInteractionEvaluatorExposer 
  : public PrecomputedInteractionEvaluator
{
public:

  static ComplexDouble ComputeElement_X(uint a, 
										uint b, 
										uint c, 
										uint d, 
										vector<vector<vector<vector<double> > > > &Vnnnn, 
										vector<vector<vector<ComplexDouble> > > & PsiB, 
										double nmax
										)
  {
	return ComputeElement(a, b, c, d, Vnnnn, PsiB, nmax);
  }

  static uint Permutations_X(uint a, 
							 uint b, 
							 uint c, 
							 uint d
							 )
  { return Permutations(a, b, c, d);}


  ComplexDouble ComputeElement_R(uint a, uint b, uint c, uint d, vector<vector<vector<vector<double> > > > &Vnnnn, vector<vector<vector<ComplexDouble> > > & PsiB, double nmax)
  {
	ComplexDouble sum = 0.0;
	for(uint n1 = 0; n1 < nmax; ++n1)
	  {
		for(uint n2 = 0; n2 < nmax; ++n2)
		  {
			for(uint n3 = 0; n3 < nmax; ++n3)
			  {
				for(uint n4 = 0; n4 < nmax; ++n4)
				  {
					sum += PsiB[0][a][n1] * 
					  PsiB[1][b][n2] * 
					  PsiB[0][c][n3] * 
					  PsiB[1][d][n4]
					  * Vnnnn[n1][n2][n3][n4];
				  }
			  }
		  }
	  }
	return sum;
  }
};







class PrecomputedInteractionEvaluatorTest : public GenericUnitTest
{
 public:
  int runUnitTests() const; ///Main function.
 protected:
  int TestCase1() const; /// Test case.
  int TestCase2() const;
};



#endif
