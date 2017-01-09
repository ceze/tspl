/*****************************************************************************
 *                               matrix_test.cpp
 *
 * Matrix class testing.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <matrix.h>


using namespace std;
using namespace splab;


typedef double  Type;
const   int     M = 3;
const   int     N = 3;


int main()
{
	Type x;

	Matrix<Type> m1;
	m1.resize( M, N );
	x = 1.0;
	m1 = x;
	cout << "matrix m1 : " << m1 << endl;

	x = 2.0;
	Matrix<Type> m2( M, N, x );
	cout << "matrix m2 : " << m2 << endl;

	Matrix<Type> m3 = m1;
	cout << "matrix m3 : " << m3 << endl;
	m3.resize( 3, 4 );
	for( int i=1; i<=3; ++i )
		for( int j=1; j<=4; ++j )
			m3(i,j) = Type(i*j);

    int row = m3.dim(1);
	int column = m3.dim(2);
	cout << "the row number of new matrix m3 : " << row << endl;
	cout << "the column number of new matrix m3 : " << column << endl;
	cout << "new matrix m3 : " << m3 << endl;

	cout << "new matrix m3 : " << m3 << endl;
	cout << "the diagonal matrix of m3 : " << diag( m3 ) << endl;
	cout << "the transpose matrix of m3 : " << trT( m3 ) << endl;

    cout << "\t\t\t\tmatrix-scalar operand" << endl << endl;
	cout << "scalar x = " << x << endl;
	cout << "m1 + x : " << m1+x << endl;
	cout << "x + m1 : " << x+m1 << endl;
	m1 += x;
	cout << "m1 += x : " << m1 << endl;
	cout << "m1 - x : " << m1-x << endl;
	cout << "x - m1 : " << x-m1 << endl;
	m1 -= x;
	cout << "m1 -= x : " << m1 << endl;
	cout << "m1 * x : " << m1*x << endl;
	cout << "x * m1 : " << x*m1 << endl;
	m1 *= x;
	cout << "m1 *= x : " << m1 << endl;
	cout << "m1 / x : " << m1/x << endl;
	cout << "x / m1 : " << x/m1 << endl;
	m1 /= x;
	cout << "m1 /= x : " << m1 << endl;

    cout << "\t\t\telementwise matrix-matrix operand" << endl << endl;
	cout << "m1 + m2 : " << m1 + m2 << endl;
	m1 += m2;
	cout << "m1 += m2 : " << m1 << endl;
	cout << "m1 - m2 : " << m1-m2 << endl;
	m1 -= m2;
	cout << "m1 -= m2 : " << m1 << endl;
	m1 *= m2;
	cout << "m1 *= m2 : " << m1 << endl;
	m1 /= m2;
	cout << "m1 /= m2 : " << m1 << endl;

	cout << "column minimum vector of m1 : " << min(m1) << endl;
	cout << "column maximum vector of m1 : " << max(m1) << endl;
	cout << "column sum vector of m1 : " << sum(m1) << endl;
	cout << "column mean vector of m1 : " << mean(m1) << endl;

    cout << "\t\t\t\tmatrix-vector operand" << endl << endl;
    Vector<Type> v1( 3, 2 );
	cout << "vector v1 : " << v1 << endl;
	cout << "m1 * v1 : " << m1*v1 << endl;
	cout << "m1^T * v1 : " << trMult(m1, v1) << endl;

	cout << "\t\t\t\tmatrix-matrix operand" << endl << endl;
	cout << "m1 * m2 : " << m1*m2 << endl;
	cout << "m3^T * m2 : " << trMult(m3, m2) << endl;
	cout << "m1 * m3^T^T : " << multTr(m1, trT(m3)) << endl;

	Matrix<int> m4( 4, 5 );
	Vector<int> v2(5);
	v2[0] = 1;  v2[1] = 2;  v2[2] = 3;  v2[3] = 4;  v2[4] = 5;
	for( int i=0; i<4; ++i )
		m4.setRow( i*v2, i );

	cout << "matrix m4 : " << m4 << endl;
	cout << "column vectors of m4 : " << endl;
	for( int j=0; j<5; ++j )
		cout << "the " << j << "th column" << m4.getColumn(j) << endl;

    v2.resize(4);
	v2[0] = 1;  v2[1] = 2;  v2[2] = 3;  v2[3] = 4;
	for( int j=0; j<5; ++j )
		m4.setColumn( j*v2, j );

	cout << "row vectors of m4 : " << endl;
	for( int i=0; i<4; ++i )
		cout << "the " << i << "th row" << m4.getRow(i) << endl;

    cout << "\t\t\t\tcomplex matrix operand" << endl << endl;
    complex<Type> c = polar(1.0,PI/4);
    Matrix< complex<Type> > A( M, N ), B(M,N), C( M, N+1 );
    for( int i=0; i<M; i++ )
		for( int j=0; j<N; j++ )
		{
		    A[i][j] = complex<Type>( Type(0.3*i+0.7*j), sin(Type(i+j)) );
		    B[i][j] = complex<Type>( cos(Type(i+j)), Type(0.8*i+0.2*j) );
		}
    for( int i=0; i<N; ++i )
        C.setColumn( A.getColumn(i), i );
    C.setColumn( B.getColumn(0), N );

    cout << setiosflags(ios::fixed) << setprecision(2);
    cout << "Matrix A:  " << A << endl;
    cout << "Matrix B:  " << B << endl;
    cout << "Matrix C:  " << C << endl;
    cout << "Absolute of C : " << abs(C) << endl;
    cout << "Angle of C : " << arg(C) << endl;
    cout << "Real part of C : " << real(C) << endl;
    cout << "Imaginary part of C : " << imag(C) << endl;

    cout << "c*A - A*B + trH(B)  " << c*A - A*B + trH(B) << endl;
    cout << "A*C * C^H*B  " << A*C * trMult(C,B) << endl;
    cout << "(A.*B+C*C^H)) ./ A^T  "
         << elemDivd( elemMult(A,B)+multTr(C,C), trT(A) ) << endl;
    cout << "diag( A * diag(C) )  " << diag( A*diag(C)/c ) * c << endl;
    cout << "A .* ( B./(c*m1) )  "
         << elemMultEq( A, elemDivdEq(B,c*complexMatrix(m1)) ) << endl;

	return 0;
}
