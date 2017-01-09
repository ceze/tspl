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
 *                            classicalpse-impl.h
 *
 * Implementation for classical power spectrum estimatoin methods.
 *
 * Zhang Ming, 2010-11, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * The correlogram power spectral estimator.
 * xn       : input signal
 * L        : the number of power spectrum density samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
inline Vector<Type> correlogramPSE( const Vector<Type> &xn, int L )
{
    return periodogramPSE( xn, rectangle(xn.size(),Type(1)), L );
}


/**
 * The windowed periodogram power spectral estimator.
 * xn       : input signal
 * wn       : window function
 * L        : the number of power spectrum density samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> periodogramPSE( const Vector<Type> &xn, const Vector<Type> &wn,
                             int L )
{
    int N = xn.size(),
	    M = wn.size();

	assert( M <= N );
	if( M < N )
	{
	    cerr << "The length of window is smaller than the length of data, ";
        cerr << "the data will be trucated to the window length!" << endl;
	}

    Vector<Type> wxn(L);
	if( L >= M )
	{
	    for( int i=0; i<M; ++i )
            wxn[i] = xn[i] * wn[i];
	}
	else
	{
        cerr << "The FFT points is smaller than the data points, ";
        cerr << "the data will be trucated to the FFT points!" << endl;
        for( int i=0; i<L; ++i )
            wxn[i] = xn[i] * wn[i];
	}

	Vector<Type> absXk = abs( fft( wxn ) );
	return absXk*absXk / Type(M);
}


/**
 * The Bartlett method of power spectral estimation.
 * xn       : input signal
 * M        : the length of subsequences
 * L        : the number of power spectrum density samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
inline Vector<Type> bartlettPSE( const Vector<Type> &xn, int M, int L )
{
    return welchPSE( xn, rectangle(M,Type(1)), M, L );
}


/**
 * The Welch method of power spectral estimation.
 * xn       : input signal
 * wn       : window function
 * K        : the number of subsequence
 * L        : the number of power spectrum density samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> welchPSE( const Vector<Type> &xn, const Vector<Type> &wn,
                       int K, int L )
{
    int N = xn.size(),
	    M = wn.size();


	assert( M < N );
	assert( K < N );

    int S = ( N-M+K )/K;
	Type P = sum( wn*wn ) / Type(M);

	Vector<Type> phi(L);
	for( int i=0; i<S; ++i )
		phi += periodogramPSE( wkeep(xn,M,i*K), wn, L );

	return phi/(S*P);
}


/**
 * The Blackman-Tukey method of power spectral estimation.
 * The correlation function is obtained from the standard biased estimate.
 * xn       : input signal
 * wn       : window function
 * L        : the number of power spectrum density samples
 * return   : spectral estimates at L frequencies:
 *            w = 0, 2*pi/L, ..., 2*pi(L-1)/L
 */
template <typename Type>
Vector<Type> btPSE( const Vector<Type> &xn, const Vector<Type> &wn, int L )
{
    int N = xn.size(),
	    M = wn.size();

	assert( M <= N );

	Vector<Type> Rxx = fastCorr( xn, "biased" );

	Vector<Type> wrn(L);
	if( L >= M )
	{
	    for( int i=0; i<M; ++i )
            wrn[i] = Rxx[N-1+i] * wn[i];
	}
	else
	{
        cerr << "The FFT points is smaller than the data points, ";
        cerr << "the data will be trucated to the FFT points!" << endl;
        for( int i=0; i<L; ++i )
            wrn[i] = Rxx[N-1+i] * wn[i];
	}

	return Type(2)*real(fft(wrn)) - wrn[0];
}
