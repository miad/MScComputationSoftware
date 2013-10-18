#ifndef SimpsonIntegrator_hpp
#define SimpsonIntegrator_hpp 1

#include <list>
using namespace std;


template<class T>
class SimpsonIntegrator
{
public:
  static T Integrate(T (*inteGrand)(double x, T k1, T k2), int nPoints, double a, double b, T k1, T k2);
private:
  SimpsonIntegrator() {} ///Thou shalt not create an instance of this class.
};


template<class T>
T SimpsonIntegrator<T>::Integrate(T (*inteGrand)(double x, T k1, T k2), int nPoints, double a, double b, T k1, T k2)
/* Integration using Simpson's rule */ 
{  
	// if n is odd - add +1 interval to make it even
	if( (nPoints/2)*2 != nPoints) 
	  ++nPoints;
    T sum = 0.0;

    double dx = (b-a)/nPoints;
    for ( int i=2; i<=nPoints-1; i=i+2 )
    {
	  double x = a+i*dx;
	  sum += 2.0*inteGrand(x,k1,k2) + 4.0*inteGrand(x+dx,k1,k2);
    }
    sum += inteGrand(a,k1,k2)+inteGrand(b,k1,k2)+4.0*inteGrand(a+dx,k1,k2);
	sum *= dx/3.0;
    return sum;
}  


#endif
