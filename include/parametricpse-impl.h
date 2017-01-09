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
 * Implementation for parametric power spectrum estimatoin methods.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * The Yule-Walker method for AR power spectral estimation.
 * xn       : input signal
 * p        : the AR model order
 * sigma2   : the variance of exciting white noise
 * return   : coefficients of AR model --- a(0), a(1), ..., a(p)
 */
template <typename Type>
Vector<Type> yulewalkerPSE( const Vector<Type> &xn, int p, Type &sigma2 )
{
    int N = xn.size();

    assert( p <= N );

    Vector<Type> rn(p+1);
    for( int i=0; i<=p; ++i )
        for( int k=0; k<N-i; ++k )
            rn[i] += xn[k+i]*xn[k];
    rn /= Type(N);

    return levinson( rn, sigma2 );
}


/**
 * The Burg method for AR power spectral estimation.
 * xn       : input signal
 * p        : the AR model order
 * sigma2   : the variance of exciting white noise
 * return   : coefficients of AR model --- a(0), a(1), ..., a(p)
 */
template <typename Type>
Vector<Type> burgPSE( const Vector<Type> &xn, int p, Type &sigma2 )
{
    int N = xn.size();
    Type numerator, denominator;
    Vector<Type> ak(p+1), akPrev(p+1), ef(N), eb(N);

    ak[0] = Type(1.0);
    sigma2 = sum(xn*xn) / Type(N);

    for( int i=1; i<N; ++i )
    {
        ef[i] = xn[i];
        eb[i-1] = xn[i-1];
    }

    for( int k=1; k<=p; ++k )
    {
        numerator = 0;
        denominator = 0;
        for( int i=k; i<N; ++i )
        {
            numerator += ef[i]*eb[i-1];
            denominator += ef[i]*ef[i] + eb[i-1]*eb[i-1];
        }
        ak[k] = -2*numerator/denominator;

        for( int i=1; i<k; ++i )
            akPrev[i] = ak[i] + ak[k]*ak[k-i];
        for( int i=1; i<k; ++i )
            ak[i] = akPrev[i];

        sigma2 *= 1 - ak[k]*ak[k];

        for( int i=N-1; i>k; --i )
        {
            ef[i] = ef[i] + ak[k]*eb[i-1];
            eb[i-1] = eb[i-2] + ak[k]*ef[i-1];
        }
    }

    return ak;
}


/**
 * Forward and backward linear prediction least square method for
 * AR power spectral estimation.
 * xn       : input signal
 * p        : the AR model order
 * sigma2   : the variance of exciting white noise
 * return   : coefficients of AR model --- a(0), a(1), ..., a(p)
 */
template <typename Type>
Vector<Type> fblplsPSE( const Vector<Type> &xn, int p, Type &sigma2 )
{
    int N = xn.size(),
        M = 2*(N-p);

    Vector<Type> u(p+1);
    u[0] = Type(1.0);

    Matrix<Type> X(M,p+1);
//    for( int i=0; i<=p; ++i )
//    {
//        for( int j=0; j<N-p; ++j )
//            X[j][i] = xn[i+j];
//        for( int j=p; j<N; ++j )
//            X[j][i] = xn[j-i];
//    }
    for( int i=0; i<N-p; ++i )
        for( int j=0; j<=p; ++j )
            X[i][j] = xn[i+j];
    for( int i=p; i<N; ++i )
        for( int j=0; j<=p; ++j )
            X[i][j] = xn[i-j];

    Matrix<Type> Rp = trMult(X,X) / Type(M);
    Vector<Type> ak = luSolver( Rp, u );
    sigma2 = 1/ak[0];
    ak *= sigma2;

    return ak;
}


/**
 * The Bartlett method of power spectral estimation.
 * ak       : AR coefficients
 * bk       : MA coefficients
 * sigma2   : the variance of exciting white noise
 * L        : the points number of PSD
 * return   : spectral density at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> armaPSD( const Vector<Type> &ak, const Vector<Type> &bk,
                      const Type &sigma2, int L )
{
    int p = ak.size()-1,
        q = bk.size()-1;
    Vector<Type> Xk(L);

    Type zRe, zIm, aRe, aIm, bRe, bIm,
         Xre, Xim,
         re, im;
	Type omega,
         den, numRe, numIm;

	for( int k=0; k<L; ++k )
	{
		omega = Type(TWOPI*k/L);
		zRe = cos(-omega);
		zIm = sin(-omega);

        // numerator
		bRe = 0;
		bIm = 0;
		for( int i=q; i>0; --i )
		{
			re = bRe;
			im = bIm;
			bRe = (re+bk[i])*zRe - im*zIm;
			bIm = (re+bk[i])*zIm + im*zRe;
		}
		bRe += bk[0];

        // denominator
		aRe = 0;
		aIm = 0;
		for( int i=p; i>0; --i )
		{
			re = aRe;
			im = aIm;
			aRe = (re+ak[i])*zRe - im*zIm;
			aIm = (re+ak[i])*zIm + im*zRe;
		}
		aRe += ak[0];

        // Power Spectrum Density
		numRe = aRe*bRe + aIm*bIm;
		numIm = aRe*bIm - aIm*bRe;
		den = aRe*aRe + aIm*aIm;
		Xre = numRe/(den);
		Xim = numIm/(den);
        Xk[k] = sigma2 * (Xre*Xre + Xim*Xim);
	}

	return Xk;
}
