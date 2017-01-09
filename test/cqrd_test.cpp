/*****************************************************************************
 *                               cqrd_test.cpp
 *
 * CQRD class testing.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <cqrd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 4;
const   int     N = 3;


int main()
{
	Matrix<Type> A(M,N);
	A[0][0] = 1;	A[0][1] = 2;	A[0][2] = 3;
	A[1][0] = 1;	A[1][1] = 2;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 9;	A[2][2] = 7;
	A[3][0] = 5;	A[3][1] = 6;	A[3][2] = 8;
//	A = trT(A);

	Matrix<complex<Type> > cA = complexMatrix( A, elemMult(A,A) );
	CQRD<Type> qr;
	qr.dec(cA);
	Matrix<complex<Type> > Q = qr.getQ();
	Matrix<complex<Type> > R = qr.getR();

    cout << setiosflags(ios::fixed) << setprecision(3);
	cout << "The original matrix cA : " << cA << endl;
	cout << "The orthogonal matrix Q  : " << Q << endl;
	cout << "The upper triangular matrix R : " << R << endl;
	cout << "Q^H * Q : " << trMult(Q,Q) << endl;
	cout << "cA - Q*R : " << cA - Q*R << endl;

    Vector<Type> b(M);
	b[0]= 1;	b[1] = 0;	b[2] = 1, b[3] = 2;
	Vector<complex<Type> > cb = complexVector(b);

	if( qr.isFullRank() )
	{
	    Vector<complex<Type> > x = qr.solve(cb);
        cout << "The constant vector cb : " << cb << endl;
        cout << "The least squares solution of cA * x = cb : " << x << endl;

        Matrix<complex<Type> > X = qr.solve( eye( M, complex<Type>(1) ) );
        cout << "The least squares solution of cA * X = I : " << X << endl;
        cout << "The cA * X: " << cA*X << endl;
	}
	else
        cout << " The matrix is rank deficient! " << endl;

	return 0;
}
