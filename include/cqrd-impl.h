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
 * Implementation for CQRD class.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
CQRD<Type>::CQRD()
{
}

template<typename Type>
CQRD<Type>::~CQRD()
{
}


/**
 * Create a QR factorization for a complex matrix A.
 */
template <typename Type>
void CQRD<Type>::dec( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols(),
        p = min(m,n);

    Type absV0, normV, beta;
    complex<Type> alpha;

    diagR = Vector<complex<Type> >(p);
    betaR = Vector<Type>(p);
	QR = A;

    // main loop.
    for( int k=0; k<p; ++k )
    {
		// Form k-th Householder vector.
		normV = 0;
		for( int i=k; i<m; ++i )
			normV += norm(QR[i][k]);
        normV = sqrt(normV);

		absV0 = abs(QR[k][k]);
		alpha = -normV * QR[k][k]/absV0;
		beta  = 1 / ( normV * (normV+absV0) );
		QR[k][k] -= alpha;

		// Apply transformation to remaining columns.
		for( int j=k+1; j<n; ++j )
		{
			complex<Type> s = 0;
			for( int i=k; i<m; ++i )
				s += conj(QR[i][k]) * QR[i][j];
            s *= beta;
			for( int i=k; i<m; ++i )
				QR[i][j] -= s*QR[i][k];
		}

		diagR[k] = alpha;
		betaR[k] = beta;
	}

//    int m = A.rows(),
//        n = A.cols(),
//        p = min(m,n);
//
//    Type absV0, normV, beta;
//    complex<Type> alpha;
//
//    diagR = Vector<complex<Type> >(p);
//    betaR = Vector<Type>(p);
//	QR = A;
//
//    // main loop.
//    for( int k=0; k<p; ++k )
//    {
//		// Form k-th Householder vector.
//		Vector<complex<Type> > v(m-k);
//		for( int i=k; i<m; ++i )
//			v[i-k] = QR[i][k];
//
//		absV0 = abs(v[0]);
//		normV = norm(v);
//		alpha = -normV * v[0]/absV0;
//		beta  = 1 / ( normV * (normV+absV0) );
//		v[0] -= alpha;
//
//		for( int i=k; i<m; ++i )
//			QR[i][k] = v[i-k];
//
//		// Apply transformation to remaining columns.
//		for( int j=k+1; j<n; ++j )
//		{
//			complex<Type> s = 0;
//			for( int i=k; i<m; ++i )
//				s += conj(v[i-k]) * QR[i][j];
//			for( int i=k; i<m; ++i )
//				QR[i][j] -= beta*s*v[i-k];
//		}
//
//		diagR[k] = alpha;
//		betaR[k] = beta;
//	}
}


/**
 * Flag to denote the matrix is of full rank.
 */
template <typename Type>
inline bool CQRD<Type>::isFullRank() const
{
    for( int j=0; j<diagR.dim(); ++j )
        if( abs(diagR[j]) == Type(0) )
            return false;

    return true;
}


/**
 * Return the upper triangular factorof the QR factorization.
 */
template <typename Type>
Matrix<complex<Type> > CQRD<Type>::getQ()
{
    int m = QR.rows(),
        p = betaR.dim();

    Matrix<complex<Type> >Q( m, p );
//	for( int i=0; i<p; ++i )
//		Q[i][i] = 1;

    for( int k=p-1; k>=0; --k )
    {
		Q[k][k] = 1 - betaR[k]*norm(QR[k][k]);
		for( int i=k+1; i<m; ++i )
			Q[i][k] = -betaR[k]*QR[i][k]*conj(QR[k][k]);

		for( int j=k+1; j<p; ++j )
		{
			complex<Type> s = 0;
			for( int i=k; i<m; ++i )
				s += conj(QR[i][k])*Q[i][j];
            s *= betaR[k];

			for( int i=k; i<m; ++i )
				Q[i][j] -= s * QR[i][k];

		}
	}

	return Q;
}


/**
 * Return the orthogonal factorof the QR factorization.
 */
template <typename Type>
Matrix<complex<Type> > CQRD<Type>::getR()
{
    int n = QR.cols(),
        p = diagR.dim();
    Matrix<complex<Type> > R( p, n );

    for( int i=0; i<p; ++i )
	{
		R[i][i] = diagR[i];

        for( int j=i+1; j<n; ++j )
            R[i][j] = QR[i][j];
	}

    return R;
}


/**
 * Least squares solution of A*x = b, where A and b are complex.
 * Return x: a n-length vector that minimizes the two norm
 * of Q*R*X-B. If B is non-conformant, or if QR.isFullRank()
 * is false, the routinereturns a null (0-length) vector.
 */
template <typename Type>
Vector<complex<Type> > CQRD<Type>::solve( const Vector<complex<Type> > &b )
{
    int m = QR.rows(),
        n = QR.cols();

    assert( b.dim() == m );

    // matrix is rank deficient
    if( !isFullRank() )
        return Vector<complex<Type> >();

    Vector<complex<Type> > x = b;

    // compute y = Q^H * b
    for( int k=0; k<n; ++k )
    {
        complex<Type> s = 0;
        for( int i=k; i<m; ++i )
            s += conj(QR[i][k])*x[i];

        s *= betaR[k];
        for( int i=k; i<m; ++i )
            x[i] -= s * QR[i][k];
    }

    // solve R*x = y;
    for( int k=n-1; k>=0; --k )
    {
        x[k] /= diagR[k];
        for( int i=0; i<k; ++i )
            x[i] -= x[k]*QR[i][k];
    }

    // return n portion of x
    Vector<complex<Type> > x_(n);
    for( int i=0; i<n; ++i )
        x_[i] = x[i];

    return x_;
}


/**
 * Least squares solution of A*X = B, where A and B are complex.
 * return X: a matrix that minimizes the two norm of Q*R*X-B.
 * If B is non-conformant, or if QR.isFullRank() is false, the
 * routinereturns a null (0) matrix.
 */
template <typename Type>
Matrix<complex<Type> > CQRD<Type>::solve( const Matrix<complex<Type> > &B )
{
    int m = QR.rows();
    int n = QR.cols();

    assert( B.rows() == m );

    // matrix is rank deficient
    if( !isFullRank() )
        return Matrix<complex<Type> >(0,0);

    int nx = B.cols();
    Matrix<complex<Type> > X = B;

    // compute Y = Q^H*B
    for( int k=0; k<n; ++k )
        for( int j=0; j<nx; ++j )
        {
            complex<Type> s = 0;
            for( int i=k; i<m; ++i )
                s += conj(QR[i][k])*X[i][j];

            s *= betaR[k];
            for( int i=k; i<m; ++i )
                X[i][j] -= s*QR[i][k];
        }

    // solve R*X = Y;
    for( int k=n-1; k>=0; --k )
    {
        for( int j=0; j<nx; ++j )
            X[k][j] /= diagR[k];

        for( int i=0; i<k; ++i )
            for( int j=0; j<nx; ++j )
                X[i][j] -= X[k][j]*QR[i][k];
    }

    // return n x nx portion of X
    Matrix<complex<Type> > X_( n, nx );
    for( int i=0; i<n; ++i )
        for( int j=0; j<nx; ++j )
            X_[i][j] = X[i][j];

     return X_;
}
