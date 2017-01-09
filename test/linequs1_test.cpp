/*****************************************************************************
 *                              linequs1_test.cpp
 *
 * Deterministic Linear Equations testing.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


#define BOUNDS_CHECK

#include <iostream>
#include <iomanip>
#include <linequs1.h>


using namespace std;
using namespace splab;


typedef float   Type;
const   int     M = 3;
const   int     N = 3;


int main()
{
	Matrix<Type> A(M,N), B;
	Vector<Type> b(N);

    // ordinary linear equations
	A[0][0] = 1;	A[0][1] = 2;	A[0][2] = 1;
	A[1][0] = 2;	A[1][1] = 5;	A[1][2] = 4;
	A[2][0] = 1;	A[2][1] = 1;	A[2][2] = 0;
	B = eye( N, Type(1.0) );
	b[0] = 1;	b[1] = 0;	b[2] = 1;
    Matrix<Type> invA( A );
    Matrix<Type> X( B );

    cout << setiosflags(ios::fixed) << setprecision(4) << endl;
    cout << "The original matrix A is : " << A << endl;
	gaussSolver( invA, X );
	cout << "The inverse of A is (Gauss Solver) : "
         << invA << endl;
    cout << "The inverse matrix of A is (LUD Solver) : "
         << luSolver( A, B ) << endl;
    cout << "The constant vector b : " << b << endl;
	cout << "The solution of A * x = b is (Gauss Solver) : "
         << gaussSolver( A, b ) << endl;
	cout << "The solution of A * x = b is (LUD Solver) : "
         << luSolver( A, b ) << endl << endl;

    Matrix<complex<Type> > cA = complexMatrix( A, B );
    Vector<complex<Type> > cb = complexVector( b, b );
    cout << "The original complex matrix cA is : " << cA << endl;
    cout << "The constant complex vector cb : " << cb << endl;
	cout << "The solution of cA * cx = cb is (Gauss Solver) : "
         << gaussSolver( cA, cb ) << endl;
	cout << "The solution of cA * cx = cb is (LUD Solver) : "
         << luSolver( cA, cb ) << endl;
    cout << "The cA*cx - cb is : "<< cA*luSolver(cA,cb) - cb << endl << endl;

    // linear equations with symmetric coefficient matrix
	for( int i=1; i<N+1; ++i )
	{
		for( int j=1; j<N+1; ++j )
			if( i == j )
				A(i,i) = Type(i);
			else
				if( i < j )
					A(i,j) = Type(i);
				else
					A(i,j) = Type(j);
		b(i) = Type( i*(i+1)/2.0 + i*(N-i) );
	}
	cout << "The original matrix A : " << A << endl;
	cout << "The inverse matrix of A is (Cholesky Solver) : "
         << choleskySolver( A, B ) << endl;
	cout << "The constant vector b : " << b << endl;
	cout << "The solution of Ax = b is (Cholesky Solver) : "
         << choleskySolver( A, b ) << endl << endl;

    cA = complexMatrix( A );
    cb = complexVector( b, b );
    cout << "The original complex matrix A : " << cA << endl;
    cout << "The constant complex vector b : " << cb << endl;
	cout << "The solution of Ax = b is (Cholesky Solver) : "
         << choleskySolver( cA, cb ) << endl << endl;

    // upper and lower triangular system
    for( int i=0; i<N; ++i )
        for( int j=i; j<N; ++j )
            B[i][j] = A[i][j];
    cout << "The original matrix B : " << B << endl;
    cout << "The constant vector b : " << b << endl;
	cout << "The solution of Ax = b is (Upper Triangular Solver) : "
         << utSolver( B, b ) << endl << endl;

    cA = complexMatrix( trT(B), -trT(B) );
    cout << "The original complex matrix A : " << cA << endl;
    cout << "The constant complex vector b : " << cb << endl;
	cout << "The solution of Ax = b is (Lower Triangular Solver) : "
         << ltSolver( cA, cb ) << endl << endl;

	// Tridiagonal linear equations
	Vector<Type> aa( 3, -1 ), bb( 4, 4 ), cc( 3, -2 ), dd( 4 );
    dd(1) = 3;   dd(2) = 2;   dd(3) = 2;   dd(4) = 3;
    cout << "The elements below main diagonal is : " << aa << endl;
    cout << "The elements on main diagonal is : " << bb << endl;
    cout << "The elements above main diagonal is : " << cc << endl;
    cout << "The elements constant vector is : " << dd << endl;
    cout << "Teh solution is (Forward Elimination and Backward Substitution) : "
         << febsSolver( aa, bb, cc, dd ) << endl;

    Vector<complex<Type> > caa = complexVector(aa);
    Vector<complex<Type> > cbb = complexVector(bb);
    Vector<complex<Type> > ccc = complexVector(cc);
    Vector<complex<Type> > cdd = complexVector(Vector<Type>(dd.dim()),dd);
    cout << "The elements below main diagonal is : " << caa << endl;
    cout << "The elements on main diagonal is : " << cbb << endl;
    cout << "The elements above main diagonal is : " << ccc << endl;
    cout << "The elements constant vector is : " << cdd << endl;
    cout << "Teh solution is (Forward Elimination and Backward Substitution) : "
         << febsSolver( caa, cbb, ccc, cdd ) << endl;

	return 0;
}
