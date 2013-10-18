#ifndef Matrix_hpp
#define Matrix_hpp

#include "basisException.hh"
#include "Globals.hpp"
#include <stdio.h>

using namespace std;

///This class should be a wrapper for a matrix stored as a 1D array. In that way, we can extract the array directly to use it with LAPACK.
template<class T>
class Matrix
{
public:
  Matrix(unsigned int _rows, ///Number of rows for the matrix to use.
		 unsigned int _columns ///Number of columns.
		 ); ///Constructor.

  ~Matrix(); ///Destructor.

  void InitializeAll(T value); ///Set all elements to

  T & Element(unsigned int row, ///Row number. Zero indexed. 
				 unsigned int column ///Column number. Zero indexed.
				 ); ///Returns the element in row row and column column. Throws exception if one or both of these parameters are out of bound.

  T * GetArray(); ///Returns the underlying array itself. Use this to send it to LAPACK. Ownership is NOT TRANSFERRED, although the array may (will) be modified (by LAPACK). Do NOT use delete on this, or you will have a segfault on your hands.

  unsigned int Rows(); ///Returns the number of rows in the matrix.

  unsigned int Columns(); ///Returns the number of columns in the matrix.

  bool IsSquare(); ///Returns true if the number of rows and the number of columns is equal.

  bool IsSymmetric(bool verbose = false ///If set to true, this will print info on the first found cells where symmetricity is violated, if any.
				   ); ///Checks if the matrix is symmetric, or not. Returns true if symmetric.
  inline bool IsHermitian(bool verbose = false ///If set to true, this will print info on the first found cells where hermiticity is violated, if any.
						  ); ///Checks if the matrix is Hermitian or not. Returns true if Hermitian.
  
private:
  T * ElementArray; ///The underlying data structure (array) containing the data.
  unsigned int rows; ///Number of rows in the array.
  unsigned int columns; ///Number of columns.
};

template<class T>
unsigned int Matrix<T>::Rows()
{
  return rows;
}

template<class T>
unsigned int Matrix<T>::Columns()
{
  return columns;
}

template<class T>
bool Matrix<T>::IsSquare()
{
  return rows == columns;
}

template<class T>
Matrix<T>::Matrix(unsigned int _rows, unsigned int _columns)
  :rows(_rows),columns(_columns)
{
  if( rows < 1 || columns < 1 )
	{
	  throw basisException("Matrix: Must have # rows > 0 and # columns > 0");
	}
  ElementArray = new T[rows*columns];
}

template<class T>
void Matrix<T>::InitializeAll(T value)
{
  for(unsigned int i = 0; i<rows*columns; ++i)
	{
	  ElementArray[i]=value;
	}
}

template<class T>
Matrix<T>::~Matrix()
{
  delete ElementArray;
}

template<class T>
bool Matrix<T>::IsSymmetric( bool verbose )
{
  if( !IsSquare() )
	{
	  return false; ///A matrix must be square in order to be symmetric.
	}
  for(unsigned int n = 0; n<rows; ++n)
	{
	  for(unsigned int m = 0; m<n; ++m)
		{
		  if(Element(n, m) != Element(m, n) )
			{
			  if(verbose)
				printf("Symmetricity invalidity detected: (%d, %d) = %lf + %lfi, (%d, %d) = %lf + %lfi.\n",n,m,real(Element(n,m)),imag(Element(n,m)),m,n,real(Element(m, n)), imag(Element(m, n)));
			}
			return false;
		}
	}
  return true;
}

template <>
inline bool Matrix<ComplexDouble>::IsHermitian(bool verbose)
{
  if ( !IsSquare() )
	{
	  return false;
	}
  for(unsigned int n = 0; n<rows; ++n)
	{
	  for(unsigned int m = 0; m<=n; ++m)
		{
		  if(Element(n, m) != conj(Element(m, n)) )
			{
			  if(verbose)
				printf("Hermiticity invalidity detected: (%d, %d) = %lf + %lfi, (%d, %d) = %lf + %lfi.\n",n,m,real(Element(n,m)),imag(Element(n,m)),m,n,real(Element(m, n)), imag(Element(m, n)));
			  return false;
			}
		}
	}
  return true;
}

template<class T>
T & Matrix<T>::Element(unsigned int row, unsigned int column)
{
  if(row >= rows)
	{
	  throw basisException("Matrix: Row out of bounds.");
	}
  if( column >= columns )
	{
	  throw basisException("Matrix Column out of bounds.");
	}
  return ElementArray[row*columns+column];
}

template<class T>
T * Matrix<T>::GetArray()
{
  return ElementArray;
}

#endif
