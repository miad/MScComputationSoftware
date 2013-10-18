#ifndef Compute_hh
#define Compute_hh 1

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <cblas.h>
#include "Potential.hh"
#include "lapacke.h"
#include "lapacke_utils.h"
#include "KPoints.hh"
#include <iostream>
#include "Globals.hpp"
#include "Matrix.hpp"
#include "SimpsonIntegrator.hpp"


using namespace std;

Potential * myPotential;


ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2);

int main(int argc, char *argv[]);

#endif
