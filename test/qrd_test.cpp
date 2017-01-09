/*****************************************************************************
 *                               qrd_test.cpp
 *
 * QRD class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <qrd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 4;
const   int     N = 3;


int main()
{
	Matrix<Type> A(M,N), Q, R;
	A[0][0] = 1;	A[0][1] = 0;	A[0][2] = 0;
	A[1][0] = 1;	A[1][1] = 2;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 3;	A[2][2] = 9;
	A[3][0] = 1;	A[3][1] = 3;	A[3][2] = 9;

	Matrix<Type> B = trT(A);
	QRD<Type> qr;
	qr.dec(B);
	Q = qr.getQ();
	R = qr.getR();

    cout << setiosflags(ios::fixed) << setprecision(4);
	cout << "The original matrix B : " << B << endl;
	cout << "The orthogonal matrix Q  : " << Q << endl;
	cout << "The upper triangular matrix R : " << R << endl;
	cout << "B - Q*R : " << B - Q*R << endl;

    Vector<Type> b(M);
	b[0]= 1;	b[1] = 0;	b[2] = 1, b[3] = 2;

	qr.dec(A);
	if( qr.isFullRank() )
	{
	    Vector<Type> x = qr.solve(b);
        cout << "The constant vector b : " << b << endl;
        cout << "The least squares solution of A * x = b : " << x << endl;

        Matrix<Type> X = qr.solve(eye(M,Type(1)));
        cout << "The least squares solution of A * X = I : " << X << endl;
        cout << "The A * X: " << A*X << endl;
	}
	else
        cout << " The matrix is rank deficient! " << endl;

	return 0;
}
