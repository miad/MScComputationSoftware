#include "Compute.hh"



ComplexDouble IntegrandValue(double x, ComplexDouble k1, ComplexDouble k2)
{
  return myPotential->Evaluate(x) * exp(ComplexDouble(0,1)*(k1-k2));
}

int main(int argc, char *argv[])
{
  int nKValues = 100;
  myPotential = new Potential();
  KPoints myKPoints(nKValues, 50, 2, 2);
  int xPoints = 2000;

  CMatrix HamiltonianMatrix(2*nKValues+1,2*nKValues+1);
  HamiltonianMatrix.InitializeAll(0.);
  
  foru(i, (int)myKPoints.GetPoints()->size())
	{
	  foru(j, (int)myKPoints.GetPoints()->size())
		{

		  HamiltonianMatrix.Element(i, j) += Integrator::Integrate( IntegrandValue, xPoints, myPotential->GetMinX(), myPotential->GetMaxX(), (*myKPoints.GetPoints())[i], (*myKPoints.GetPoints())[j] );

		  if ( i==j )
			HamiltonianMatrix.Element(i, j) += pow(HBAR, 2)/(2.*MASS) * pow((*myKPoints.GetPoints())[i], 2);
		}
	}
  delete myPotential;
  if ( ! HamiltonianMatrix.IsHermitian(true) )
	{
	  throw basisException("The matrix was found to be non-hermitian.");
	}

  /*  InitializePotential();
  Complex * test = new Complex[9];
  foru(i, 9)
	test[i]=new Complex(i, i);

  //GetXPoints(NUMBER_OF_X_SAMPLE_POINTS, XMIN, XMAX);
  //GetKPoints(NUMBER_OF_K_VALUES, K_CUTOFF, TRIANGLE_KMID, TRIANGLE_KDEPTH);

  printf("Matrix:\n");

  foru(i, 3)
	{
	  printf("(");
	  foru(j, 3)
		{
		  printf("%lf %lf   ", real(test[3*i+j]), imag(test[3*i+j]));
		}
	  printf(")\n");
	}


  int matrix_order = LAPACK_ROW_MAJOR;
  char jobvl = 'N';
  char jobvr = 'N';
  lapack_int n = 3;
  cplex * a = test;
  lapack_int lda = 3;
  lapack_complex_double * w = new Complex[n];
  lapack_complex_double * vl = 0;
  lapack_int ldvl = 3;
  lapack_complex_double * vr = 0;
  lapack_int ldvr = 3;

  int reply = LAPACKE_zgeev(matrix_order, jobvl, jobvr,
						   n, a, lda, w, vl, ldvl, vr, ldvr);

  printf("Reply value: %d\n", reply);


  
  foru(i, 3)
	  printf("Eigenvalue: %lf + %lfi\n", real(w[i]), imag(w[i]));


lapack_int LAPACKE_zgeev( int matrix_order, char jobvl, char jobvr,
                          lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* w,
                          lapack_complex_double* vl, lapack_int ldvl,
                          lapack_complex_double* vr, lapack_int ldvr );

							
  ClearPotential();
*/
  return 0;
}
