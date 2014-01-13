#ifndef HermiteRule_hh
#define HermiteRule_hh 1

# include <cstdlib>
# include <cmath>
# include <ctime>
# include <cstring>
#include <iostream>
#include <vector>
#include <utility>
#include "RLMacros.hpp"
#include "RLException.hh"

using namespace std;

class HermiteRule
{
public:
  static vector<pair<double, double> > GetRule(unsigned int numberOfPoints, ///Number of points in the rule.
											   double a = -1, ///Center point for rule.
											   double b = 1 ///Scaling for rule.
											   ); ///Computes the Hermite rule for n points: 
  //Integral ( -oo < x < +oo ) f(x) exp ( - b * ( x - a )^2 ) dx

private:
  HermiteRule() { } ///Prevent construction of object.


  ///Not checked.
  static void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
			   double wts[] );
  static void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
			  double t[], double wts[] );
  static double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
						double bj[] );
  static void imtqlx ( int n, double d[], double e[], double z[] );
  static void parchk ( int kind, int m, double alpha, double beta );
  static double r8_abs ( double x );
  static double r8_epsilon ( );
  static double r8_gamma ( double x );
  static double r8_huge ( );
  static double r8_sign ( double x );
  static void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
					 double swts[], double st[], int kind, double alpha, double beta, double a, 
					 double b );
  static void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
			  double wts[] );
};

#endif
