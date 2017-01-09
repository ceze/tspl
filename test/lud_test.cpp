/*****************************************************************************
 *                               lud_test.cpp
 *
 * LUD class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <lud.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 4;


int main()
{
	Matrix<Type> A(3,3), B(3,3), L, U;
	Vector<Type> b(3);

	A[0][0] = 1;	A[0][1] = 2;	A[0][2] = 1;
	A[1][0] = 2;	A[1][1] = 5;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 1;	A[2][2] = 0;

	B[0][0] = 1;	B[0][1] = 0;	B[0][2] = 0;
	B[1][0] = 0;	B[1][1] = 1;	B[1][2] = 0;
	B[2][0] = 0;	B[2][1] = 0;	B[2][2] = 1;

	b[0] = 1;	b[1] = 0;	b[2] = 1;

    LUD<Type> lu;
	lu.dec(A);
	L = lu.getL();
	U = lu.getU();
    cout << setiosflags(ios::fixed) << setprecision(4);
	cout << "The original matrix A is : " << A << endl;
	cout << "The unit lower triangular matrix L is : " << L << endl;
	cout << "The upper triangular matrix U is : " << U << endl;

	Vector<Type> x = lu.solve(b);
	cout << "The constant vector b : " << b << endl;
	cout << "The solution of A * x = b : " << x << endl;
	cout << "The solution of A * x - b : " << A*x - b << endl;

	cout << "The inverse matrix of A : " << lu.solve(B) << endl;
	cout << "The A * inverse(A) : " << A*lu.solve(B) << endl;
	cout << "The determinant of A : " << endl;
	cout << lu.det() << endl << endl << endl;

    Matrix<complex<Type> > cA(M,N), cPA(M,N), invcA, cL, cU;
    for( int i=0; i<M; i++ )
		for( int j=0; j<N; j++ )
			cA[i][j] = complex<Type>( Type(0.3*i+0.7*j), sin(Type(i+j)) );
	LUD<complex<Type> > clu;

	clu.dec(cA);
	cL = clu.getL();
	cU = clu.getU();
	Vector<int> p = clu.getPivot();
	for( int i=0; i<M; ++i )
        cPA.setRow( cA.getRow(p[i]), i );

    cout << setiosflags(ios::fixed) << setprecision(3);
	cout << "The original complex matrix cA is : " << cA << endl;
	cout << "The unit lower triangular matrix cL is : " << cL << endl;
	cout << "The upper triangular matrix cU is : " << cU << endl;
	cout << "cP*cA - cL*cU is : " << cPA - cL*cU << endl;

	if( M == N )
	{
	    invcA = clu.solve( eye( M, complex<Type>(1.0,0) ) );
	    cout << "The inverse matrix of cA : " << invcA << endl;
	    cout << "The cA * inverse(cA) : " << cA*invcA << endl;
	}

	return 0;
}
