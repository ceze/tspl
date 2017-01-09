/*****************************************************************************
 *                               csvd_test.cpp
 *
 * CSVD class testing.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <csvd.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 4;
const   int     N = 5;


int main()
{
    Matrix< complex<Type> > A( M, N ), B(M,N);
    for( int i=0; i<M; i++ )
		for( int j=0; j<N; j++ )
			A[i][j] = complex<Type>( Type(0.3*i+0.7*j), sin(Type(i+j)) );

	CSVD<Type> svd;
	svd.dec(A);

	Matrix< complex<Type> > U = svd.getU();
	Matrix< complex<Type> > V = svd.getV();
	Matrix<Type> S = svd.getSM();
    Matrix< complex<Type> > CS = complexMatrix( S );
//	Vector<Type> S = svd.getSV();

    cout << setiosflags(ios::fixed) << setprecision(3);
	cout << "Matrix--------A: " << A << endl;
	cout << "Matrix--------U: " << U << endl;
	cout << "Vector--------S: " << S << endl;
	cout << "Matrix--------V: " << V << endl;
    cout << "Matrix--------A - U * S * V^H:  "
         << A - U*CS*trH(V) << endl;
//         << A - U*multTr(CS,V) << endl;

	cout << "The rank of A : " << svd.rank() << endl << endl;
	cout << "The condition number of A : " << svd.cond() << endl << endl;

	return 0;
}
