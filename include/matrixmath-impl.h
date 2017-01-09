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
 *                               vector-impl.h
 *
 * Implementation for Matrix math functions.
 *
 * Zhang Ming, 2010-08, Xi'an Jiaotong University.
 *****************************************************************************/


template <typename Type>
Matrix<Type> abs( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = abs( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> cos( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = cos( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> sin( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = sin( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> tan( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = tan( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> acos( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = acos( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> asin( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = asin( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> atan( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = atan( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> exp( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = exp( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> log( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = log( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> log10( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = log10( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> sqrt( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = sqrt( A[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> pow( const Matrix<Type> &B, const Matrix<Type> &E )
{
    int m = B.rows(),
        n = B.cols();
    assert( m == E.rows() );
    assert( n == E.cols() );

    Matrix<Type> tmp( m, n );
    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = pow( B[i][j], E[i][j] );

    return tmp;
}


template <typename Type>
Matrix<Type> pow( const Matrix<Type> &B, const Type &e )
{
    int m = B.rows(),
        n = B.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = pow( B[i][j], e );

    return tmp;
}


template <typename Type>
Matrix<Type> pow( const Type &b, const Matrix<Type> &E )
{
    int m = E.rows(),
        n = E.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = pow( b, E[i][j] );

    return tmp;
}
