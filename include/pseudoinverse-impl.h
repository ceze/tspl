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
 *                             pseudoinverse-impl.h
 *
 * Implementation for matrix pseudoinverse
 *
 * Zhang Ming, 2010-08 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Compute the pseudoinverse of a real matrix.
 */
template <typename Real>
Matrix<Real> pinv( const Matrix<Real> &A, Real tol )
{
    int m = A.rows(),
        n = A.cols();

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();

    int r = 0;
    if( tol <= 0 )
        tol = max( m, n ) * s[0] * EPS;
    for( int i=0; i<s.size(); ++i )
        if( s[i] >= tol )
            r++;

    for( int i=0; i<n; ++i )
        for( int k=0; k<r; ++k )
                V[i][k] /= s[k];

    Matrix<Real> invA( n, m );
    for( int i=0; i<n; ++i )
        for( int j=0; j<m; ++j )
        {
            Real sum = 0;
            for( int k=0; k<r; ++k )
                sum += V[i][k]*U[j][k];
            invA[i][j] = sum;
        }

    return invA;
}


/**
 * Compute the pseudoinverse of a complex matrix.
 */
template <typename Type>
Matrix<complex<Type> > pinv( const Matrix<complex<Type> > &A, Type tol )
{
    int m = A.rows(),
        n = A.cols();

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();

    int r = 0;
    if( tol <= 0 )
        tol = max( m, n ) * s[0] * EPS;
    for( int i=0; i<s.size(); ++i )
        if( s[i] >= tol )
            r++;

    for( int i=0; i<n; ++i )
        for( int k=0; k<r; ++k )
                V[i][k] /= s[k];

    Matrix<complex<Type> >invA( n, m );
    for( int i=0; i<n; ++i )
        for( int j=0; j<m; ++j )
        {
            complex<Type> sum = 0;
            for( int k=0; k<r; ++k )
                sum += V[i][k]*conj(U[j][k]);
            invA[i][j] = sum;
        }

    return invA;
}
