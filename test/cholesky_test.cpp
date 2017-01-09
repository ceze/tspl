/*****************************************************************************
 *                               cholesky_test.cpp
 *
 * Cholesky class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <cholesky.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 5;


int main()
{
	Matrix<Type> A(N,N), L(N,N);
	Vector<Type> b(N);

	for( int i=1; i<N+1; ++i )
	{
		for( int j=1; j<N+1; ++j )
			if( i == j )
				A(i,i) = i;
			else
				if( i < j )
					A(i,j) = i;
				else
					A(i,j) = j;

		b(i) = i*(i+1)/2.0 + i*(N-i);
	}

	cout << setiosflags(ios::fixed) << setprecision(3);
	cout << "The original matrix A : " << A << endl;
	Cholesky<Type> cho;
    cho.dec(A);

	if( !cho.isSpd() )
		cout << "Factorization was not complete." << endl;
	else
    {
        L = cho.getL();
        cout << "The lower triangular matrix L is : " << L << endl;
        cout << "A - L*L^T is : " << A - L*trT(L) << endl;

        Vector<Type> x = cho.solve(b);
        cout << "The constant vector b : " << b << endl;
        cout << "The solution of Ax = b : " << x << endl;
        cout << "The Ax - b : " << A*x-b << endl;

        Matrix<Type> IA = cho.solve(eye(N,Type(1)));
        cout << "The inverse matrix of A : " << IA << endl;
        cout << "The product of  A*inv(A) : " << A*IA << endl;
    }

	return 0;
}
