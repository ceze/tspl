/*****************************************************************************
 *                              inverse_test.cpp
 *
 * Matrix inverse testing.
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <inverse.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 3;

int main()
{
	Matrix<Type> A, invA, B(N,N);
	A.resize(3,3);
	A[0][0] = 1;	A[0][1] = 2;	A[0][2] = 1;
	A[1][0] = 2;	A[1][1] = 5;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 1;	A[2][2] = 0;

	cout << setiosflags(ios::fixed) << setprecision(4);
	cout << "The original matrix A is : " << A << endl;
	invA = inv(A);
	cout << "The inverse matrix of A (LUD) : " << invA << endl;
    invA = colPivInv(A);
    cout << "The inverse matrix of A (column pivot) : " << invA << endl;
    invA = cmpPivInv(A);
    cout << "The invese matrix of A (complete pivot) : " << invA << endl;
    cout << "The multiplication of A and its inverse: " << A*invA << endl;

    for( int i=1; i<=N; ++i )
	{
		for( int j=1; j<=N; ++j )
			if( i == j )
				B(i,i) = i;
			else if( i < j )
                B(i,j) = i;
            else
                B(i,j) = j;
	}
	cout << "The original matrix B is : " << B << endl;
	invA = inv(B,"spd");
    cout << "The inverse matrix of B (Cholesky) : " << invA << endl;
    invA = colPivInv(B);
    cout << "The inverse matrix of B (column pivot) : " << invA << endl;
    invA = cmpPivInv(B);
    cout << "The inverse matrix of B (complete pivot) : " << invA << endl;
    cout << "The multiplication of B and its inverse: " << B*invA << endl;

    cout << setiosflags(ios::fixed) << setprecision(3);
    Matrix<complex<Type> > cA(N,N), cIA;
    cA = complexMatrix( A, B );

    cIA = inv(cA);
    cout << "The original complex matrix cA is: " << cA << endl;
    cout << "The inverse matrix of cA is (general method): "
         << inv(cA) << endl;
    cout << "The inverse matrix of cA is (column pivot): "
         << colPivInv(cA) << endl;
    cout << "The inverse matrix of cA is (complete pivot): "
         << cmpPivInv(cA) << endl;
    cout << "The inverse matrix of cA is (real inverse): "
         << cinv(cA) << endl;
    cout << "The multiplication of cA and its inverse: "
         << cA*cIA << endl;

	return 0;
}
