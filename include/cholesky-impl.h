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
 *                               cholesky-impl.h
 *
 * Implementation for Cholesky class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
Cholesky<Type>::Cholesky() : spd(true)
{
}

template<typename Type>
Cholesky<Type>::~Cholesky()
{
}


/**
 * return true, if original matrix is symmetric positive-definite.
 */
template<typename Type>
inline bool Cholesky<Type>::isSpd() const
{
    return spd;
}


/**
 * Constructs a lower triangular matrix L, such that L*L'= A.
 * If A is not symmetric positive-definite (SPD), only a partial
 * factorization is performed. If isspd() evalutate true then
 * the factorizaiton was successful.
 */
template <typename Type>
void Cholesky<Type>::dec( const Matrix<Type> &A )
{
    int m = A.rows();
    int n = A.cols();

    spd = (m == n);
    if( !spd )
        return;

    L = Matrix<Type>(n,n);

    // main loop.
    for( int j=0; j<n; ++j )
    {
        Type d = 0;
        for( int k=0; k<j; ++k )
        {
            Type s = 0;
            for( int i=0; i<k; ++i )
                s += L[k][i]*L[j][i];

            L[j][k] = s = (A[j][k]-s) / L[k][k];
            d = d + s*s;
            spd = spd && (A[k][j] == A[j][k]);
        }

        d = A[j][j] - d;
        spd = spd && ( d > 0 );

        L[j][j] = sqrt( d > 0 ? d : 0 );
        for( int k=j+1; k<n; ++k )
            L[j][k] = 0;
    }
}


/**
 * return the lower triangular factor, L, such that L*L'=A.
 */
template<typename Type>
inline Matrix<Type> Cholesky<Type>::getL() const
{
    return L;
}


/**
 * Solve a linear system A*x = b, using the previously computed
 * cholesky factorization of A: L*L'.
 */
template <typename Type>
Vector<Type> Cholesky<Type>::solve( const Vector<Type> &b )
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

    // solve L'*x = y
    for( int k=n-1; k>=0; --k )
    {
        for( int i=k+1; i<n; ++i )
            x[k] -= x[i]*L[i][k];

        x[k] /= L[k][k];
    }

    return x;
}


/**
 * Solve a linear system A*X = B, using the previously computed
 * cholesky factorization of A: L*L'.
 */
template <typename Type>
Matrix<Type> Cholesky<Type>::solve( const Matrix<Type> &B )
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

    // solve L'*X = Y
    for( int j=0; j<nx; ++j )
        for( int k=n-1; k>=0; --k )
        {
            for( int i=k+1; i<n; ++i )
                X[k][j] -= X[i][j]*L[i][k];

            X[k][j] /= L[k][k];
        }

    return X;
}


/**
 * Main loop of specialized member function. This macro definition is
 * aimed at avoiding code duplication.
 */
#define MAINLOOP                                        \
{                                                       \
    int m = A.rows();                                   \
    int n = A.cols();                                   \
                                                        \
    spd = (m == n);                                     \
    if( !spd )                                          \
        return;                                         \
                                                        \
    for( int j=0; j<A.rows(); ++j )                     \
    {                                                   \
        spd = spd && (imag(A[j][j]) == 0);              \
        d = 0;                                          \
                                                        \
        for( int k=0; k<j; ++k )                        \
        {                                               \
            s = 0;                                      \
            for( int i=0; i<k; ++i )                    \
                s += L[k][i] * conj(L[j][i]);           \
                                                        \
            L[j][k] = s = (A[j][k]-s) / L[k][k];        \
            d = d + norm(s);                            \
            spd = spd && (A[k][j] == conj(A[j][k]));    \
        }                                               \
                                                        \
        d = real(A[j][j]) - d;                          \
        spd = spd && ( d > 0 );                         \
                                                        \
        L[j][j] = sqrt( d > 0 ? d : 0 );                \
        for( int k=j+1; k<A.rows(); ++k )               \
            L[j][k] = 0;                                \
    }                                                   \
}


/**
 * Solving process of specialized member function. This macro definition is
 * aimed at avoiding code duplication.
 */
#define SOLVE1                                          \
{                                                       \
    for( int k=0; k<n; ++k )                            \
    {                                                   \
        for( int i=0; i<k; ++i )                        \
            x[k] -= x[i]*L[k][i];                       \
                                                        \
        x[k] /= L[k][k];                                \
    }                                                   \
                                                        \
    for( int k=n-1; k>=0; --k )                         \
    {                                                   \
        for( int i=k+1; i<n; ++i )                      \
            x[k] -= x[i]*conj(L[i][k]);                 \
                                                        \
        x[k] /= L[k][k];                                \
    }                                                   \
                                                        \
    return x;                                           \
}


/**
 * Solving process of specialized member function. This macro definition is
 * aimed at avoiding code duplication.
 */
#define SOLVE2                                          \
{                                                       \
    int nx = B.cols();                                  \
    for( int j=0; j<nx; ++j )                           \
        for( int k=0; k<n; ++k )                        \
        {                                               \
            for( int i=0; i<k; ++i )                    \
                X[k][j] -= X[i][j]*L[k][i];             \
                                                        \
            X[k][j] /= L[k][k];                         \
        }                                               \
                                                        \
    for( int j=0; j<nx; ++j )                           \
        for( int k=n-1; k>=0; --k )                     \
        {                                               \
            for( int i=k+1; i<n; ++i )                  \
                X[k][j] -= X[i][j]*conj(L[i][k]);       \
                                                        \
            X[k][j] /= L[k][k];                         \
        }                                               \
                                                        \
    return X;                                           \
}


/**
 * Specializing for "dec" member function.
 */
template <>
void Cholesky<complex<float> >::dec( const Matrix<complex<float> > &A )
{
    float d;
    complex<float> s;
    L = Matrix<complex<float> >(A.cols(),A.cols());

    MAINLOOP;
}

template <>
void Cholesky<complex<double> >::dec( const Matrix<complex<double> > &A )
{
    double d;
    complex<double> s;
    L = Matrix<complex<double> >(A.cols(),A.cols());

    MAINLOOP;
}

template <>
void Cholesky<complex<long double> >::dec( const Matrix<complex<long double> > &A )
{
    long double d;
    complex<long double> s;
    L = Matrix<complex<long double> >(A.cols(),A.cols());

    MAINLOOP;
}


/**
 * Specializing for "solve" member function.
 */
template <>
Vector<complex<float> > Cholesky<complex<float> >::solve( const Vector<complex<float> > &b )
{
    int n = L.rows();
    if( b.dim() != n )
        return Vector<complex<float> >();
    Vector<complex<float> > x = b;

    SOLVE1
}

template <>
Vector<complex<double> > Cholesky<complex<double> >::solve( const Vector<complex<double> > &b )
{
    int n = L.rows();
    if( b.dim() != n )
        return Vector<complex<double> >();
    Vector<complex<double> > x = b;

    SOLVE1
}

template <>
Vector<complex<long double> > Cholesky<complex<long double> >::solve( const Vector<complex<long double> > &b )
{
    int n = L.rows();
    if( b.dim() != n )
        return Vector<complex<long double> >();
    Vector<complex<long double> > x = b;

    SOLVE1
}


/**
 * Specializing for "solve" member function.
 */
template <>
Matrix<complex<float> > Cholesky<complex<float> >::solve( const Matrix<complex<float> > &B )
{
    int n = L.rows();
    if( B.rows() != n )
        return Matrix<complex<float> >();
    Matrix<complex<float> > X = B;

    SOLVE2
}

template <>
Matrix<complex<double> > Cholesky<complex<double> >::solve( const Matrix<complex<double> > &B )
{
    int n = L.rows();
    if( B.rows() != n )
        return Matrix<complex<double> >();
    Matrix<complex<double> > X = B;

    SOLVE2
}

template <>
Matrix<complex<long double> > Cholesky<complex<long double> >::solve( const Matrix<complex<long double> > &B )
{
    int n = L.rows();
    if( B.rows() != n )
        return Matrix<complex<long double> >();
    Matrix<complex<long double> > X = B;

    SOLVE2
}
