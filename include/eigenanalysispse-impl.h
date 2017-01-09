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
 *                            parametricpse-impl.h
 *
 * Implementation for eigenanalysis algorithms for spectrum estimation.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * The Capon method (minimum variance method) for spectral estimation.
 * xn       : input signal
 * M        : the order of the covariance matrix
 * L        : the number of estimated spectral samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> caponPSE( const Vector<Type> &xn, int M, int L )
{
    int N = xn.size();

    assert( M < N );

    Vector<Type> rn(M), Ere(M), Eim(M), Tre(M), Tim(M),
                 tn(M), wn(M), Px(L);

    // auto-correlation matrix R
    for( int i=0; i<M; ++i )
    {
        tn[i] = Type(i);
        for( int k=0; k<N-i; ++k )
            rn[i] += xn[k+i]*xn[k];
    }
    rn /= Type(N-M);

    for( int k=0; k<L; ++k )
    {
        // real part and imaginary part of E(omega)
        wn = Type( TWOPI*k/L ) * tn;
        Ere = cos( wn );
        Eim = sin( wn );

        // Tre = inv(R)*Ere, Tim = inv(R)*Eim. Because R is Toeplitz, this can
        // be solved through Levinson algorithm.
        Tre = levinson( rn, Ere );
        Tim = levinson( rn, Eim );

        // denominator of spectrum
        Px[k] = dotProd(Ere,Tre) + dotProd(Eim,Tim);
    }

    return Type(M)/Px;
}


/**
 * The MUSIC method for spectral estimation.
 * xn       : input signal
 * M        : the order of the covariance matrix
 * p        : the model order
 * L        : the number of estimated spectral samples
 * return   : spectral estimates (dB) at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> musicPSE( const Vector<Type> &xn, int M, int p, int L )
{
    int N = xn.size();

    assert( M < N );
    assert( p < M );

    Type  Tre, Tim, sum;
    Vector<Type> rm(M), tm(M),
                 Ere(M), Eim(M), Vi(M),
                 wk(M), Px(L);

    // auto-correlation matrix R
    for( int i=0; i<M; ++i )
    {
        tm[i] = Type(i);
        for( int k=0; k<N-i; ++k )
            rm[i] += xn[k+i]*xn[k];
    }
    rm /= Type(N-M);
    Matrix<Type> Rx = toeplitz(rm);

    SVD<Type> svd;
    svd.dec(Rx);
    Matrix<Type> U = svd.getU();

    for( int k=0; k<L; ++k )
    {
        // real part and imaginary part of E(omega)
        wk = Type( TWOPI*k/L ) * tm;
        Ere = cos( wk );
        Eim = sin( wk );

        // Tre = Ere*Vi, Tim = Eim*Vi;
        sum = 0;
        for( int i=p; i<M; ++i )
        {
            Vi = U.getColumn(i);
            Tre = dotProd( Ere, Vi );
            Tim = dotProd( Eim, Vi );
            sum += Tre*Tre + Tim*Tim;
        }

        // spectrum
//        Px[k] = 1 / ( sum );
        Px[k] = -10*log10(sum);
    }

    return Px;
}


/**
 * The Pisarenko method for spectral estimation.
 * xn       : input signal
 * M        : the order of the covariance matrix
 * p        : the model order
 * L        : the number of estimated spectral samples
 * return   : spectral estimates (dB) at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> pisarenkoPSE( const Vector<Type> &xn, int M, int p, int L )
{
    int N = xn.size();

    assert( M < N );
    assert( p < M );

    Type  Tre, Tim;
    Vector<Type> rm(M), tm(M),
                 Ere(M), Eim(M), Vi(M),
                 wk(M), Px(L);

    // auto-correlation matrix R
    for( int i=0; i<M; ++i )
    {
        tm[i] = Type(i);
        for( int k=0; k<N-i; ++k )
            rm[i] += xn[k+i]*xn[k];
    }
    rm /= Type(N-M);
    Matrix<Type> Rx = toeplitz(rm);

    SVD<Type> svd;
    svd.dec(Rx);
    Matrix<Type> U = svd.getU();

    for( int k=0; k<L; ++k )
    {
        // real part and imaginary part of E(omega)
        wk = Type( TWOPI*k/L ) * tm;
        Ere = cos( wk );
        Eim = sin( wk );

        // Tre = Ere*Vi, Tim = Eim*Vi;
        Vi = U.getColumn(p);
        Tre = dotProd( Ere, Vi );
        Tim = dotProd( Eim, Vi );

        // spectrum
//        Px[k] = 1 / ( Tre*Tre + Tim*Tim );
        Px[k] = -10*log10(Tre*Tre+Tim*Tim);
    }

    return Px;
}


/**
 * The ESPRIT method for spectral estimation.
 * xn       : input signal
 * M        : the order of the covariance matrix
 * p        : the model order
 * return   : the esitmated frequency in the interval [-0.5,0.5]
 */
template <typename Type>
Vector<Type> espritPSE( Vector<Type> &xn, int M, int p )
{
    int N = xn.size();

    assert( M < N );
    assert( p < M );

    Vector<Type> rm(M), fk(p);
    Matrix<Type> S1(M-1,p), S2(M-1,p);

    // get the auto-correlation matrix
    for( int i=0; i<M; ++i )
        for( int k=0; k<N-i; ++k )
            rm[i] += xn[k+i]*xn[k];
    rm /= Type(N-M);
    Matrix<Type> Rx = toeplitz(rm);

    // get the eigendecomposition of R, use svd because it sorts eigenvalues
    SVD<Type> svd;
    svd.dec(Rx);
    Matrix<Type> U = svd.getU();

    // compute S1 and S2
    for( int i=0; i<M-1; ++i )
        for( int j=0; j<p; ++j )
        {
            S1[i][j] = U[i][j];
            S2[i][j] = U[i+1][j];
        }

    // compute matrix Phi
    Matrix<Type> Phi = choleskySolver( trMult(S1,S1), trMult(S1,S2) );

    // compute eigenvalues of Phi
    EVD<Type> evd;
    evd.dec(Phi);
    Vector<Type> evRe = real( evd.getCD() );
    Vector<Type> evIM = imag( evd.getCD() );

    // compute normalized frequency in the interval [-0.5,0.5]
    for( int i=0; i<p; ++i )
        fk[i] = atan2( evIM[i], evRe[i] ) / Type(TWOPI);

    return fk;
}


/**
 * The model order estimation based on minimizing the MDL criterion.
 * xn       : input signal
 * M        : the order of the covariance matrix
 * return   : p ---> the number of sinusoids
 */
template <typename Type>
int orderEst( const Vector<Type> &xn, int M )
{
    int N = xn.size();

    assert( M < N );

    int p = 0;
    Type Gp, Ap, Ep, MDL, minMDL;
    Vector<Type> rm(M);

    // auto-correlation matrix R
    for( int i=0; i<M; ++i )
        for( int k=0; k<N-i; ++k )
            rm[i] += xn[k+i]*xn[k];
    rm /= Type(N-M);
    Matrix<Type> Rx = toeplitz(rm);

    SVD<Type> svd;
    svd.dec(Rx);
    Vector<Type> S = svd.getSV();

    minMDL = Type( 0.5*(M-1)*(M+1)*log10(1.0*N) );
    for( int i=0; i<M-2; ++i )
    {
        // compute MDL(i)
        Gp = Type(1);
        Ap = Type(0);
        Ep = Type( 0.5*i*(2*M-i)*log10(1.0*N) );
        for( int j=i+1; j<M; ++j )
        {
            Gp *= S[j];
            Ap += S[j];
        }
        Ap = pow( Ap/(M-i), Type(M-i) );
        MDL = -N*log10(Gp/Ap) + Ep;

        // find the minimum MDL(i)
        if( MDL < minMDL )
        {
            p = i;
            minMDL = MDL;
        }
    }

    return p;
}
