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
 *                               qrd-impl.h
 *
 * Implementation for QRD class.
 *
 * Zhang Ming, 2010-01 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Real>
QRD<Real>::QRD()
{
}

template<typename Real>
QRD<Real>::~QRD()
{
}


/**
 * Create a QR factorization object for A.
 */
template <typename Real>
void QRD<Real>::dec( const Matrix<Real> &A )
{

    int m = A.rows(),
        n = A.cols(),
        p = min(m,n);
    QR = A;
    RDiag = Vector<Real>(p);

    // main loop.
    for( int k=0; k<p; ++k )
    {
        // Compute 2-norm of k-th column without under/overflow.
        Real nrm = 0;
        for( int i=k; i<m; ++i )
            nrm = hypot( nrm, QR[i][k] );
//        for( int i=k; i<m; ++i )
//            nrm += QR[i][k]*QR[i][k];
//        nrm = sqrt(nrm);

        if( nrm != 0 )
        {
            // Form k-th Householder vector.
            if( QR[k][k] < 0 )
                nrm = -nrm;

            for( int i=k; i<m; ++i )
                QR[i][k] /= nrm;

            QR[k][k] += 1;

            // Apply transformation to remaining columns.
            for( int j=k+1; j<n; ++j )
            {
                Real s = 0;
                for( int i=k; i<m; ++i )
                    s += QR[i][k]*QR[i][j];

                s = -s/QR[k][k];
                for( int i=k; i<m; ++i )
                    QR[i][j] += s*QR[i][k];
            }
        }

        RDiag[k] = -nrm;
    }
}


/**
 * Flag to denote the matrix is of full rank.
 */
template <typename Real>
inline bool QRD<Real>::isFullRank() const
{
    for( int j=0; j<RDiag.dim(); ++j )
        if( RDiag[j] == 0 )
            return false;

    return true;
}


/**
 * Return the upper triangular factorof the QR factorization.
 */
template <typename Real>
Matrix<Real> QRD<Real>::getQ()
{
    int m = QR.rows(),
        p = RDiag.dim();
    Matrix<Real> Q( m, p );

    for( int k=p-1; k>=0; --k )
    {
        for( int i=0; i<m; ++i )
            Q[i][k] = 0;

        Q[k][k] = 1;
        for( int j=k; j<p; ++j )
            if( QR[k][k] != 0 )
            {
                Real s = 0;
                for( int i=k; i<m; ++i )
                    s += QR[i][k] * Q[i][j];

                s = -s / QR[k][k];
                for( int i=k; i<m; ++i )
                    Q[i][j] += s*QR[i][k];
            }
    }

    return Q;
}


/**
 * Return the orthogonal factorof the QR factorization.
 */
template <typename Real>
Matrix<Real> QRD<Real>::getR()
{
    int n = QR.cols(),
        p = RDiag.dim();
    Matrix<Real> R( p, n );

    for( int i=0; i<p; ++i )
        for( int j=0; j<n; ++j )
            if( i < j )
                R[i][j] = QR[i][j];
            else if( i == j )
                R[i][j] = RDiag[i];

    return R;
}


/**
 * Retreive the Householder vectors from QR factorization
 */
template <typename Real>
Matrix<Real> QRD<Real>::getH()
{
    int m = QR.rows(),
        p = RDiag.dim();
    Matrix<Real> H( m, p );

    for( int i=0; i<m; ++i )
        for( int j=0; j<=i&&j<p; ++j )
            H[i][j] = QR[i][j];

    return H;
}


/**
 * Least squares solution of A*x = b
 * Return x: a vector that minimizes the two norm of Q*R*X-B.
 * If B is non-conformant, or if QR.isFullRank() is false,
 * the routinereturns a null (0-length) vector.
 */
template <typename Real>
Vector<Real> QRD<Real>::solve( const Vector<Real> &b )
{
    int m = QR.rows(),
        n = QR.cols();

    assert( b.dim() == m );

    // matrix is rank deficient
    if( !isFullRank() )
        return Vector<Real>();

    Vector<Real> x = b;

    // compute y = transpose(Q)*b
    for( int k=0; k<n; ++k )
    {
        Real s = 0;
        for( int i=k; i<m; ++i )
            s += QR[i][k]*x[i];

        s = -s/QR[k][k];
        for( int i=k; i<m; ++i )
            x[i] += s*QR[i][k];
    }

    // solve R*x = y;
    for( int k=n-1; k>=0; --k )
    {
        x[k] /= RDiag[k];
        for( int i=0; i<k; ++i )
            x[i] -= x[k]*QR[i][k];
    }

    // return n portion of x
    Vector<Real> x_(n);
    for( int i=0; i<n; ++i )
        x_[i] = x[i];

    return x_;
}


/**
 * Least squares solution of A*X = B
 * return X: a matrix that minimizes the two norm of Q*R*X-B.
 * If B is non-conformant, or if QR.isFullRank() is false, the
 * routinereturns a null (0) matrix.
 */
template <typename Real>
Matrix<Real> QRD<Real>::solve( const Matrix<Real> &B )
{
    int m = QR.rows();
    int n = QR.cols();

    assert( B.rows() == m );

    // matrix is rank deficient
    if( !isFullRank() )
        return Matrix<Real>(0,0);

    int nx = B.cols();
    Matrix<Real> X = B;

    // compute Y = transpose(Q)*B
    for( int k=0; k<n; ++k )
        for( int j=0; j<nx; ++j )
        {
            Real s = 0;
            for( int i=k; i<m; ++i )
                s += QR[i][k]*X[i][j];

            s = -s/QR[k][k];
            for( int i=k; i<m; ++i )
                X[i][j] += s*QR[i][k];
        }

    // solve R*X = Y;
    for( int k=n-1; k>=0; --k )
    {
        for( int j=0; j<nx; ++j )
            X[k][j] /= RDiag[k];

        for( int i=0; i<k; ++i )
            for( int j=0; j<nx; ++j )
                X[i][j] -= X[k][j]*QR[i][k];
    }

    // return n x nx portion of X
    Matrix<Real> X_( n, nx );
    for( int i=0; i<n; ++i )
        for( int j=0; j<nx; ++j )
            X_[i][j] = X[i][j];

     return X_;
}
