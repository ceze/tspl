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
 *                               linequs2-impl.h
 *
 * Implementation for solving overdetermined and underdetermined linear
 * equations.
 *
 * Zhang Ming, 2010-07 (revised 2010-12), Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Overdetermined linear equationequations solution by Least Squares
 * Generalized Inverse.
 * A  --->  The m-by-n(m>n) coefficient matrix(Full Column Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 least squares solution vector.
 */
template <typename Type>
Vector<Type> lsSolver( const Matrix<Type> &A, const Vector<Type> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() > A.cols() );

    Cholesky<Type> cho;
    cho.dec( trMult(A,A) );
    if( cho.isSpd() )
        return cho.solve( trMult(A,b) );
    else
        return luSolver( trMult(A,A), trMult(A,b) );
}


/**
 * Overdetermined linear equationequations solution by QR Decomposition.
 * A  --->  The m-by-n(m>n) coefficient matrix(Full Column Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 least squares solution vector.
 */
template <typename Real>
Vector<Real> qrLsSolver( const Matrix<Real> &A, const Vector<Real> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() > A.cols() );

    QRD<Real> qr;
    qr.dec( A );
    if( !qr.isFullRank() )
    {
        cerr << "The matrix A is not Full Rank!" << endl;
        return Vector<Real>();
    }
    else
        return qr.solve( b );
}


/**
 * Overdetermined linear equationequations solution by SVD Decomposition.
 * A  --->  The m-by-n(m>n) coefficient matrix(Full Column Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 least squares solution vector.
 */
template<typename Real>
Vector<Real> svdLsSolver( const Matrix<Real> &A, const Vector<Real> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() > A.cols() );

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();

//    for( int i=0; i<V.rows(); ++i )
//        for( int k=0; k<s.dim(); ++k )
//            V[i][k] /= s[k];
//
//    return V * trMult(U,b);

    return V * ( trMult(U,b) / s );
}


/**
 * Overdetermined linear equationequations solution by QR Decomposition.
 * A  --->  The m-by-n(m>n) coefficient matrix(Full Column Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 least squares solution vector.
 */
template<typename Type>
Vector<complex<Type> > qrLsSolver( const Matrix<complex<Type> > &A,
                                   const Vector<complex<Type> > &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() > A.cols() );

    CQRD<Type> qr;
    qr.dec( A );
    if( !qr.isFullRank() )
    {
        cerr << "The matrix A is not Full Rank!" << endl;
        return Vector<complex<Type> >();
    }
    else
        return qr.solve( b );
}


/**
 * Overdetermined linear equationequations solution by SVD Decomposition.
 * A  --->  The m-by-n(m>n) coefficient matrix(Full Column Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 least squares solution vector.
 */
template<typename Type>
Vector<complex<Type> > svdLsSolver( const Matrix<complex<Type> > &A,
                                    const Vector<complex<Type> > &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() > A.cols() );

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();

//    for( int i=0; i<V.rows(); ++i )
//        for( int k=0; k<s.dim(); ++k )
//            V[i][k] /= s[k];
//
//    return V * trMult(U,b);

    return V * ( trMult(U,b) / complexVector(s) );
}


/**
 * Undetermined linear equationequations solution by Minimum Norm
 * Generalized Inverse.
 * A  --->  The m-by-n(m<n) coefficient matrix(Full Row Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 minimum norm solution vector.
 */
template <typename Type>
Vector<Type> lnSolver( const Matrix<Type> &A, const Vector<Type> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() < A.cols() );

    Cholesky<Type> cho;
    cho.dec( multTr(A,A) );
    if( cho.isSpd() )
        return trMult( A, cho.solve(b) );
    else
        return trMult( A, luSolver(multTr(A,A),b) );
}


/**
 * Undetermined linear equationequations solution by QR Decomposition.
 * A  --->  The m-by-n(m<n) coefficient matrix(Full Row Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 minimum norm solution vector.
 */
template <typename Real>
Vector<Real> qrLnSolver( const Matrix<Real> &A, const Vector<Real> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() < A.cols() );

    Matrix<Real> At( trT( A ) );
    QRD<Real> qr;
    qr.dec( At );
    if( !qr.isFullRank() )
    {
        cerr << "The matrix A is not Full Rank!" << endl;
        return Vector<Real>();
    }
    else
    {
        Matrix<Real> Q, R;
        Q = qr.getQ();
        R = qr.getR();
        Vector<Real> y( ltSolver( trT( R ), b ) );
        return Q * y;
    }
}


/**
 * Undetermined complex linear equationequations solution by SVD
 * Decomposition.
 * A  --->  The m-by-n(m<n) coefficient matrix(Full Row Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 minimum norm solution vector.
 */
template<typename Real>
Vector<Real> svdLnSolver( const Matrix<Real> &A, const Vector<Real> &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() < A.cols() );

    SVD<Real> svd;
    svd.dec( A );
    Matrix<Real> U = svd.getU();
    Matrix<Real> V = svd.getV();
    Vector<Real> s = svd.getSV();

//    for( int i=0; i<V.rows(); ++i )
//        for( int k=0; k<s.dim(); ++k )
//            V[i][k] /= s[k];
//
//    return V * trMult(U,b);

    return V * ( trMult(U,b) / s );
}


/**
 * Undetermined complex linear equationequations solution by QR
 * Decomposition.
 * A  --->  The m-by-n(m<n) coefficient matrix(Full Row Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 minimum norm solution vector.
 */
template<typename Type>
Vector<complex<Type> > qrLnSolver( const Matrix<complex<Type> > &A,
                                   const Vector<complex<Type> > &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() < A.cols() );

    Matrix<complex<Type> > At( trH( A ) );
    CQRD<Type> qr;
    qr.dec( At );
    if( !qr.isFullRank() )
    {
        cerr << "The matrix A is not Full Rank!" << endl;
        return Vector<complex<Type> >();
    }
    else
    {
        Matrix<complex<Type> > Q, R;
        Q = qr.getQ();
        R = qr.getR();
        Vector<complex<Type> > y( ltSolver( trH( R ), b ) );
        return Q * y;
    }
}


/**
 * Undetermined complex linear equationequations solution by SVD
 * Decomposition.
 * A  --->  The m-by-n(m<n) coefficient matrix(Full Row Rank);
 * b  --->  The n-by-1 right-hand side vector;
 * x  --->  The n-by-1 minimum norm solution vector.
 */
template<typename Type>
Vector<complex<Type> > svdLnSolver( const Matrix<complex<Type> > &A,
                                    const Vector<complex<Type> > &b )
{
    assert( A.rows() == b.size() );
    assert( A.rows() < A.cols() );

    CSVD<Type> svd;
    svd.dec( A );
    Matrix<complex<Type> > U = svd.getU();
    Matrix<complex<Type> > V = svd.getV();
    Vector<Type> s = svd.getSV();

//    for( int i=0; i<V.rows(); ++i )
//        for( int k=0; k<s.dim(); ++k )
//            V[i][k] /= s[k];
//
//    return V * trMult(U,b);

    return V * ( trMult(U,b) / complexVector(s) );
}
