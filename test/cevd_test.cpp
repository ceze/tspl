/*****************************************************************************
 *                               cevd_test.cpp
 *
 * CEVD class testing.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <cevd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     N = 4;


int main()
{
	Matrix<Type> B(N,N);
    B[0][0] = 3.0;	B[0][1] = -2.0;	B[0][2] = -0.9;     B[0][3] = 0.0;
	B[1][0] = -2.0;	B[1][1] = 4.0;	B[1][2] = 1.0;      B[1][3] = 0.0;
	B[2][0] = 0.0;	B[2][1] = 0.0;	B[2][2] = -1.0;     B[2][3] = 0.0;
	B[3][0] = -0.5;	B[3][1] = -0.5;	B[3][2] = 0.1;      B[3][3] = 1.0;

	Matrix<complex<Type> > A = complexMatrix( B, elemMult(B,B) );
//	A = multTr(A,A);
	cout << setiosflags(ios::fixed) << setprecision(2);
    cout << "The original complex matrix A : " << A << endl;

	CEVD<Type> eig;
	eig.dec(A);

    if( eig.isHertimian() )
    {
        Matrix<complex<Type> > V = eig.getV();
        Vector<Type> D = eig.getRD();
        Matrix<complex<Type> > DM = diag( complexVector(D) );
        cout << "The eigenvectors matrix V is : " << V << endl;
        cout << "The eigenvalue D is : " << diag(D) << endl;
        cout << "The V'*V : " << trMult(V,V) << endl;
        cout << "The A*V - V*D : " << A*V - V*DM << endl;
    }
    else
    {
        Matrix<complex<Type> > V = eig.getV();
        Vector<complex<Type> > D = eig.getD();
        Matrix<complex<Type> > DM = diag( D );
        cout << "The complex eigenvectors matrix V : " << V << endl;
        cout << "The complex eigenvalue D : " << DM << endl;
        cout << "The A*V - V*D : " << A*V - V*DM << endl;
    }

	return 0;
}
