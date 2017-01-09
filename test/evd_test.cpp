/*****************************************************************************
 *                               evd_test.cpp
 *
 * EVD class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <evd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 2;


int main()
{
	Matrix<Type> B(N,N);
    B[0][0] = 3.0;	B[0][1] = -2.0;
	B[1][0] = -2.0;	B[1][1] = 4.0;

    Matrix<Type> A(2*N,2*N), D, V;
    Matrix<complex<Type> > cA(2*N,2*N), cD, cV;
    for( int i=0; i<N; ++i )
        for( int j=0; j<N; ++j )
        {
            A[i][j]        = B[i][j];
            A[i][j+N]      = -B[j][i];
            A[i+N][j]      = B[j][i];
            A[i+N][j+N]    = B[i][j];
        }
//    A = B;
    cA = complexMatrix(A);
    cout << setiosflags(ios::fixed) << setprecision(2);
    cout << "The original matrix A : " << A << endl;

	EVD<Type> eig;
	eig.dec(A);
	if( !eig.isComplex() )
	{
	    V = eig.getV();
	    D = diag(eig.getD());
        cout << "The eigenvectors matrix V : " << V << endl;
        cout << "The eigenvalue D : " << D << endl;
        cout << "The V'*V : " << trMult(V,V) << endl;
        cout << "The A*V - V*D : " << A*V - V*D << endl;
	}
	else
	{
	    cV = eig.getCV();
	    cD = diag(eig.getCD());
        cout << "The complex eigenvectors matrix V : " << cV << endl;
        cout << "The complex eigenvalue D : " << cD << endl;
        cout << "The A*V - V*D : " << cA*cV - cV*cD << endl;
	}

	return 0;
}
