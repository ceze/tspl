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
 *                               cevd-impl.h
 *
 * Implementation for CEVD class.
 *
 * Zhang Ming, 2010-12, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructor and destructor
 */
template<typename Type>
CEVD<Type>::CEVD() : hermitian(true)
{
}

template<typename Type>
CEVD<Type>::~CEVD()
{
}


/**
 * Check for symmetry, then construct the eigenvalue decomposition
 */
template <typename Type>
void CEVD<Type>::dec( const Matrix<complex<Type> > &A )
{
    int N = A.cols();

    assert( A.rows() == N );

    V = Matrix<complex<Type> >(N,N);

    Matrix<Type> S(2*N,2*N);
    for( int i=0; i<N; ++i )
        for( int j=0; j<N; ++j )
        {
            S[i][j]     = A[i][j].real();
            S[i][j+N]   = -A[i][j].imag();
            S[i+N][j]   = -S[i][j+N];
            S[i+N][j+N] = S[i][j];
        }

    EVD<Type> eig;
    eig.dec(S);

    if( eig.isSymmetric() )
    {
        hermitian = true;
        rd = Vector<Type>(N);

        Matrix<Type> RV = eig.getV();
        Vector<Type> Rd = eig.getD();

        for( int i=0; i<N; ++i )
            rd[i] = Rd[2*i];

        for( int j=0; j<N; ++j )
        {
            int j2 = 2*j;
            for( int i=0; i<N; ++i )
                V[i][j] = complex<Type>( RV[i][j2], RV[i+N][j2] );
        }
    }
    else
    {
        hermitian = false;
        d = Vector<complex<Type> >(N);

        Matrix<complex<Type> > cV = eig.getCV();
        Vector<complex<Type> > cd = eig.getCD();

        for( int i=0; i<N; ++i )
            d[i] = cd[2*i];

        for( int j=0; j<N; ++j )
        {
            int j2 = 2*j;
            for( int i=0; i<N; ++i )
                V[i][j] = complex<Type>( cV[i][j2].real()-cV[i+N][j2].imag(),
                                         cV[i][j2].imag()+cV[i+N][j2].real() );
        }
    }
}


/**
 * If the matrix is Hermitian, then return true.
 */
template <typename Type>
bool CEVD<Type>::isHertimian() const
{
    return hermitian;
}


/**
 * Return the COMPLEX eigenvector matrix
 */
template <typename Type>
inline Matrix<complex<Type> > CEVD<Type>::getV() const
{
    return V;
}


/**
 * Return the complex eigenvalues vector.
 */
template <typename Type>
inline Vector<complex<Type> > CEVD<Type>::getD() const
{
    return d;
}


/**
 * Return the real eigenvalues vector.
 */
template <typename Type>
inline Vector<Type> CEVD<Type>::getRD() const
{
    return rd;
}
