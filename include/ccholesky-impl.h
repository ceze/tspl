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
 *                             ccholesky-impl.h
 *
 * Implementation for CCholesky (complex Cholesky) class.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
CCholesky<Type>::CCholesky() : spd(true)
{
}

template<typename Type>
CCholesky<Type>::~CCholesky()
{
}


/**
 * return true, if original matrix is symmetric positive-definite.
 */
template<typename Type>
inline bool CCholesky<Type>::isSpd() const
{
    return spd;
}


/**
 * Constructs a lower triangular matrix L, such that L*L^H= A.
 * If A is not symmetric positive-definite (SPD), only a partial
 * factorization is performed. If isspd() evalutate true then
 * the factorizaiton was successful.
 */
template <typename Type>
void CCholesky<Type>::dec( const Matrix<Type> &A )
{
    int m = A.rows();
    int n = A.cols();

    spd = (m == n);
    if( !spd )
        return;

    L = Matrix<Type>(A.cols(),A.cols());

    // main loop
    for( int j=0; j<A.rows(); ++j )
    {
        Type d = 0;
        spd = spd && (imag(A[j][j]) == 0);

        for( int k=0; k<j; ++k )
        {
            Type s = 0;
            for( int i=0; i<k; ++i )
                s += L[k][i] * conj(L[j][i]);

            L[j][k] = s = (A[j][k]-s) / L[k][k];
            d = d + s*conj(s);
            spd = spd && (A[k][j] == conj(A[j][k]));
        }

        d = A[j][j] - d;
        spd = spd && ( real(d) > 0 );

        L[j][j] = sqrt( real(d) > 0 ? d : 0 );
        for( int k=j+1; k<A.rows(); ++k )
            L[j][k] = 0;
    }
}


/**
 * return the lower triangular factor, L, such that L*L'=A.
 */
template<typename Type>
inline Matrix<Type> CCholesky<Type>::getL() const
{
    return L;
}


/**
 * Solve a linear system A*x = b, using the previously computed
 * cholesky factorization of A: L*L'.
 */
template <typename Type>
Vector<Type> CCholesky<Type>::solve( const Vector<Type> &b )
{
    int n = L.rows();
    if( b.dim() != n )
        return Vector<Type>();

    Vector<Type> x = b;

    // solve L*y = b
    for( int k=0; k<n; ++k )
    {
        for( int i=0; i<k; ++i )
            x[k] -= x[i]*L[k][i];

        x[k] /= L[k][k];
    }

    // solve L^H*x = y
    for( int k=n-1; k>=0; --k )
    {
        for( int i=k+1; i<n; ++i )
            x[k] -= x[i]*conj(L[i][k]);

        x[k] /= L[k][k];
    }

    return x;
}


/**
 * Solve a linear system A*X = B, using the previously computed
 * cholesky factorization of A: L*L'.
 */
template <typename Type>
Matrix<Type> CCholesky<Type>::solve( const Matrix<Type> &B )
{
    int n = L.rows();
    if( B.rows() != n )
        return Matrix<Type>();

    Matrix<Type> X = B;
    int nx = B.cols();

    // solve L*Y = B
    for( int j=0; j<nx; ++j )
        for( int k=0; k<n; ++k )
        {
            for( int i=0; i<k; ++i )
                X[k][j] -= X[i][j]*L[k][i];

            X[k][j] /= L[k][k];
        }

    // solve L^H*x = y
    for( int j=0; j<nx; ++j )
        for( int k=n-1; k>=0; --k )
        {
            for( int i=k+1; i<n; ++i )
                X[k][j] -= X[i][j]*conj(L[i][k]);

            X[k][j] /= L[k][k];
        }

    return X;
}
