/*****************************************************************************
 *                              linequs2_test.cpp
 *
 * Undetermined Linear Equations testing.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <linequs2.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 3;


int main()
{
	Matrix<Type> A(M,N), B(M,N);
	Vector<Type> b(N);

    // overdetermined linear equations
    A.resize( 4, 3 );
	A[0][0] = 1;	A[0][1] = -1;	A[0][2] = 1;
	A[1][0] = 1;	A[1][1] = 2;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 3;	A[2][2] = 9;
	A[3][0] = 1;	A[3][1] = -4;	A[3][2] = 16;
	b.resize( 4 );
	b[0]= 1;    b[1] = 2;   b[2] = 3;   b[3] = 4;

	cout << setiosflags(ios::fixed) << setprecision(3);
    cout << "The original matrix A : " << A << endl;
    cout << "The constant vector b : " << b << endl;
	cout << "The least square solution is (using generalized inverse) : "
         << lsSolver( A, b ) << endl;
    cout << "The least square solution is (using QR decomposition) : "
         << qrLsSolver( A, b ) << endl;
    cout << "The least square solution is (using SVD decomposition) : "
         << svdLsSolver( A, b ) << endl;

    Matrix<complex<Type> > cA = complexMatrix( A, A );
    Vector<complex<Type> > cb = complexVector( b );
    cout << "The original complex matrix cA : " << cA << endl;
    cout << "The constant complex vector cb : " << cb << endl;
	cout << "The least square solution is (using generalized inverse) : "
         << lsSolver( cA, cb ) << endl;
    cout << "The least square solution is (using QR decomposition) : "
         << qrLsSolver( cA, cb ) << endl;
    cout << "The least square solution is (using SVD decomposition) : "
         << svdLsSolver( cA, cb ) << endl;

    // undetermined linear equations
    Matrix<Type> At( trT( A ) );
	b.resize( 3 );
	b[0]= 1;	b[1] = 2;   b[2]= 3;
	cout << "The original matrix A : " << At << endl;
    cout << "The constant vector b : " << b << endl;
	cout << "The least norm solution is (using generalized inverse) : "
         << lnSolver( At, b ) << endl;
	cout << "The least norm solution is (using QR decomposition) : "
         << qrLnSolver( At, b ) << endl;
    cout << "The least norm solution is (using SVD decomposition) : "
         << svdLnSolver( At, b ) << endl;

    cA = complexMatrix( At, -At );
    cb = complexVector( b );
    cout << "The original complex matrix cA : " << cA << endl;
    cout << "The constant complex vector cb : " << cb << endl;
	cout << "The least square solution is (using generalized inverse) : "
         << lnSolver( cA, cb ) << endl;
    cout << "The least square solution is (using QR decomposition) : "
         << qrLnSolver( cA, cb ) << endl;
    cout << "The least square solution is (using SVD decomposition) : "
         << svdLnSolver( cA, cb ) << endl;

	return 0;
}
