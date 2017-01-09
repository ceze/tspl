/*
 * Copyright (c) 2008-2011 Zhang Ming (M. Zhang), zmjerry@163.com
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 2 or any later version.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details. A copy of the GNU General Public License is available at:
 * http://www.fsf.org/licensing/licenses
 */


/*****************************************************************************
 *                            linequs1-impl.h
 *
 * Implementation for deterministic linear equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Linear equations solution by Gauss-Jordan elimination, Ax=B.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<Type> gaussSolver( const Matrix<Type> &A, const Vector<Type> &b )
{
    assert( A.rows() == b.size() );

    int rows = b.size();
    Vector<Type> x( rows );
    Matrix<Type> C(A);
    Matrix<Type> B( rows, 1 );
    for( int i=0; i<rows; ++i )
        B[i][0] = b[i];

    gaussSolver( C, B );

    for( int i=0; i<rows; ++i )
        x[i] = B[i][0];

    return x;
}


/**
 * Linear equations solution by Gauss-Jordan elimination, AX=B.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * B  --->  The n-by-m right-hand side vectors;
 * When solving is completion, A is replaced by its inverse and
 * B is replaced by the corresponding set of solution vectors.
 * Adapted from Numerical Recipes.
 */
template <typename Type>
void gaussSolver( Matrix<Type> &A, Matrix<Type> &B )
{
    assert( A.rows() == B.rows() );

    int i, j, k, l, ll, icol, irow,
        n=A.rows(),
        m=B.cols();
    Type big, dum, pivinv;

    Vector<int> indxc(n), indxr(n), ipiv(n);
    for( j=0; j<n; j++ )
        ipiv[j]=0;

    for( i=0; i<n; ++i )
    {
        big = 0.0;
        for( j=0; j<n; ++j )
            if( ipiv[j] != 1 )
                for( k=0; k<n; ++k )
                    if( ipiv[k] == 0 )
                        if( abs(A[j][k]) >= abs(big) )
                        {
                            big = abs(A[j][k]);
                            irow = j;
                            icol = k;
                        }

        ++(ipiv[icol]);

        if( irow != icol )
        {
            for( l=0; l<n; ++l )
                swap( A[irow][l], A[icol][l] );
            for( l=0; l<m; ++l )
                swap( B[irow][l], B[icol][l] );
        }
        indxr[i] = irow;
        indxc[i] = icol;

        if( abs(A[icol][icol]) == 0.0 )
        {
            cerr << "Singular Matrix!" << endl;
            return;
        }

        pivinv = Type(1) / A[icol][icol];
        A[icol][icol] = Type(1);
        for( l=0; l<n; ++l )
            A[icol][l] *= pivinv;
        for( l=0; l<m; ++l )
            B[icol][l] *= pivinv;
        for( ll=0; ll<n; ++ll )
            if( ll != icol )
            {
                dum = A[ll][icol];
                A[ll][icol] = 0;
                for( l=0; l<n; ++l )
                    A[ll][l] -= A[icol][l]*dum;
                for( l=0; l<m; ++l )
                    B[ll][l] -= B[icol][l]*dum;
            }
    }

    for( l=n-1; l>=0; --l )
        if( indxr[l] != indxc[l] )
            for( k=0; k<n; k++ )
                swap( A[k][indxr[l]], A[k][indxc[l]] );
}


/**
 * Linear equations solution by LU decomposition, Ax=b.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<Type> luSolver( const Matrix<Type> &A, const Vector<Type> &b )
{
    assert( A.rows() == b.size() );

    LUD<Type> lu;
    lu.dec( A );
    return lu.solve( b );
}


/**
 * Linear equations solution by LU decomposition, AX=B.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * B  --->  The n-by-m right-hand side vector;
 * X  --->  The n-by-m solution vectors.
 */
template <typename Type>
Matrix<Type> luSolver( const Matrix<Type> &A, const Matrix<Type> &B )
{
    assert( A.rows() == B.rows() );

    LUD<Type> lu;
    lu.dec( A );
    return lu.solve( B );
}


/**
 * Linear equations solution by Cholesky decomposition, Ax=b.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<Type> choleskySolver( const Matrix<Type> &A, const Vector<Type> &b )
{
    assert( A.rows() == b.size() );

    Cholesky<Type> cho;
    cho.dec(A);
    if( cho.isSpd() )
        return cho.solve( b );
    else
    {
        cerr << "Factorization was not complete!" << endl;
        return Vector<Type>(0);
    }
}


/**
 * Linear equations solution by Cholesky decomposition, AX=B.
 * A  --->  The n-by-n coefficient matrix(Full Rank);
 * B  --->  The n-by-m right-hand side vector;
 * X  --->  The n-by-m solution vectors.
 */
template <typename Type>
Matrix<Type> choleskySolver( const Matrix<Type> &A, const Matrix<Type> &B )
{
    assert( A.rows() == B.rows() );

    Cholesky<Type> cho;
    cho.dec(A);
    if( cho.isSpd() )
        return cho.solve( B );
    else
    {
        cerr << "Factorization was not complete!" << endl;
        return Matrix<Type>(0,0);
    }
}


/**
 * Solve the upper triangular system U*x = b.
 * U  --->  The n-by-n upper triangular matrix(Full Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
*/
template <typename Type>
Vector<Type> utSolver( const Matrix<Type> &U, const Vector<Type> &b )
{
    int n = b.dim();

    assert( U.rows() == n );
    assert( U.rows() == U.cols() );

	Vector<Type> x(b);
	for( int k=n; k >= 1; --k )
	{
		x(k) /= U(k,k);
		for( int i=1; i<k; ++i )
			x(i) -= x(k) * U(i,k);
    }

    return x;
}


/**
 * Solve the lower triangular system L*x = b.
 * L  --->  The n-by-n lower triangular matrix(Full Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
*/
template <typename Type>
Vector<Type> ltSolver( const Matrix<Type> &L, const Vector<Type> &b )
{
	int n = b.dim();

    assert( L.rows() == n );
    assert( L.rows() == L.cols() );

	Vector<Type> x(b);
	for( int k=1; k <= n; ++k )
	{
		x(k) /= L(k,k);
		for( int i=k+1; i<= n; ++i )
			x(i) -= x(k) * L(i,k);
    }

    return x;
}


/**
 * Tridiagonal equations solution by Forward Elimination and Backward Substitution.
 * a  --->  The n-1-by-1 main(0th) diagonal vector;
 * b  --->  The n-by-1   -1th(above main) diagonal vector;
 * c  --->  The n-1-by-1 +1th(below main) diagonal vector;
 * d  --->  The right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<Type> febsSolver( const Vector<Type> &a, const Vector<Type> &b,
                         const Vector<Type> &c, const Vector<Type> &d )
{
    int n = b.size();

    assert( a.size() == n-1 );
    assert( c.size() == n-1 );
    assert( d.size() == n );

    Type mu = 0;
    Vector<Type> x(d),
                 bb(b);

    for( int i=0; i<n-1; ++i )
    {
        mu = a[i] / bb[i];
        bb[i+1] = bb[i+1] - mu*c[i];
        x[i+1] = x[i+1] - mu*x[i];
    }

    x[n-1] = x[n-1] / bb[n-1];
    for( int i=n-2; i>=0; --i )
        x[i] = ( x[i]-c[i]*x[i+1] ) / bb[i];

//        for( int j=1; j<n; ++j )
//        {
//            mu = a(j) / bb(j);
//            bb(j+1) = bb(j+1) - mu*c(j);
//            x(j+1) = x(j+1) - mu*x(j);
//        }
//
//        x(n) = x(n) / bb(n);
//        for( int j=n-1; j>0; --j )
//            x(j) = ( x(j)-c(j)*x(j+1) ) / bb(j);

    return x;
}
