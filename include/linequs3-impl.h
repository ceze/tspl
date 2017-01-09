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
 *                            linequs3-impl.h
 *
 * Implementation for Rank Defect linear equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Rank defect linear equationequations solution by Truncated SVD.
 * A  --->  The m-by-n coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Real>
Vector<Real> tsvd( const Matrix<Real> &A, const Vector<Real> &b, Real tol )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();
    Vector<Real> x(n);

    int r = 0;
    if( tol <= 0 )
        tol = max( m, n ) * s[0] * EPS;
    for( int i=0; i<s.size(); ++i )
        if( s[i] >= tol )
            r++;

    // y = U^T * b
    Vector<Real> y(r);
    for( int i=0; i<r; ++i )
        for( int j=0; j<m; ++j )
            y[i] += U[j][i]*b[j];

    // y = y / s
    for( int i=0; i<r; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<r; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<r; ++j )
//        {
//            Real sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += U[k][j]*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}


/**
 * Rank defect complex linear equationequations solution by Truncated SVD.
 * A  --->  The m-by-n coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<complex<Type> > tsvd( const Matrix<complex<Type> > &A,
                             const Vector<complex<Type> > &b,
                             Type tol )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();
    Vector<complex<Type> > x(n);

    int r = 0;
    if( tol <= 0 )
        tol = max( m, n ) * s[0] * EPS;
    for( int i=0; i<s.size(); ++i )
        if( s[i] >= tol )
            r++;

    // y = U^H * b
    Vector<complex<Type> > y(r);
    for( int i=0; i<r; ++i )
        for( int j=0; j<m; ++j )
            y[i] += conj(U[j][i])*b[j];

    // y = y / s
    for( int i=0; i<r; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<r; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<r; ++j )
//        {
//            complex<Type> sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += conj(U[k][j])*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}


/**
 * Rank defect linear equationequations solution by Dampted SVD.
 * A  --->  The m-by-n(m>n) coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * sigma :  dampted factor;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Real>
Vector<Real> dsvd( const Matrix<Real> &A, const Vector<Real> &b,
                   Real &sigma )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();
    Vector<Real> x(n);

    int p = s.size();
    s += sigma;

    // y = U^T * b
    Vector<Real> y(p);
    for( int i=0; i<p; ++i )
        for( int j=0; j<m; ++j )
            y[i] += U[j][i]*b[j];

    // y = y / s
    for( int i=0; i<p; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<p; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<p; ++j )
//        {
//            Real sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += U[k][j]*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}


/**
 * Rank defect complex linear equationequations solution by Dampted SVD.
 * A  --->  The m-by-n(m>n) coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * sigma :  dampted factor;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<complex<Type> > dsvd( const Matrix<complex<Type> > &A,
                             const Vector<complex<Type> > &b,
                             Type &sigma )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();
    Vector<complex<Type> > x(n);

    int p = s.size();
    s += sigma;

    // y = U^H * b
    Vector<complex<Type> > y(p);
    for( int i=0; i<p; ++i )
        for( int j=0; j<m; ++j )
            y[i] += conj(U[j][i])*b[j];

    // y = y / s
    for( int i=0; i<p; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<p; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<p; ++j )
//        {
//            complex<Type> sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += conj(U[k][j])*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}


/**
 * Rank defect linear equationequations solution by Tikhonov Regularization.
 * A  --->  The m-by-n(m>n) coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * alpha :  regularization factor;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Real>
Vector<Real> tikhonov( const Matrix<Real> &A, const Vector<Real> &b,
                       Real &alpha )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();
    Vector<Real> x(n);

    int p = s.size();
    Real alpha2 = alpha*alpha;
    for( int i=0; i<p; ++i )
        s[i] += alpha2/s[i];

    // y = U^T * b
    Vector<Real> y(p);
    for( int i=0; i<p; ++i )
        for( int j=0; j<m; ++j )
            y[i] += U[j][i]*b[j];

    // y = y / s
    for( int i=0; i<p; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<p; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<p; ++j )
//        {
//            Real sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += U[k][j]*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}


/**
 * Rank defect complex linear equationequations solution by
 * Tikhonov Regularization.
 * A  --->  The m-by-n(m>n) coefficient matrix(Rank Defect);
 * b  --->  The n-by-1 right-hand side vector;
 * alpha :  regularization factor;
 * x  --->  The n-by-1 solution vector.
 */
template <typename Type>
Vector<complex<Type> > tikhonov( const Matrix<complex<Type> > &A,
                                 const Vector<complex<Type> > &b,
                                 Type &alpha )
{
    int m = A.rows(),
        n = A.cols();

    assert( m == b.size() );

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();
    Vector<complex<Type> > x(n);

    int p = s.size();
    Type alpha2 = alpha*alpha;
    for( int i=0; i<p; ++i )
        s[i] += alpha2/s[i];

    // y = U^H * b
    Vector<complex<Type> > y(p);
    for( int i=0; i<p; ++i )
        for( int j=0; j<m; ++j )
            y[i] += conj(U[j][i])*b[j];

    // y = y / s
    for( int i=0; i<p; ++i )
        y[i] /= s[i];

    // x = V * y
    for( int i=0; i<n; ++i )
        for( int j=0; j<p; ++j )
            x[i] += V[i][j]*y[j];

//    for( int i=0; i<n; ++i )
//        for( int j=0; j<p; ++j )
//        {
//            complex<Type> sum = 0;
//            for( int k=0; k<m; ++k )
//                sum += conj(U[k][j])*b[k];
//            x[i] += sum * V[i][j] / s[j];
//        }

    return x;
}
