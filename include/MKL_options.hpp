#ifndef MKL_OPTIONS_HPP
#define MKL_OPTIONS_HPP 1



#ifdef USE_MKL_LAPACKE
#include <complex>
#include "mkl_lapacke.h"
#define RIKARD_COMPLEX_TYPE MKL_Complex16
#define RIKARD_COMPLEX_COPY_TO(from, to){				\
	to.real = from.real();								\
	to.imag = from.imag();								\
  }
#define RIKARD_COMPLEX_COPY_FROM(from, to){				\
	to.real() = from.real;								\
	to.imag() = from.imag;								\
  }
#else
#define RIKARD_COMPLEX_TYPE ComplexDouble
#include "lapacke.h"
#include "lapacke_utils.h"
#define RIKARD_COMPLEX_COPY_TO(from, to){		\
	to = from;									\
  }
#define RIKARD_COMPLEX_COPY_FROM(from, to){		\
	to = from;									\
  }
#endif


#endif
