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
 *                                wvd-impl.h
 *
 * Implementation for Wigner-Ville distribution.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Compute WVD of 1D real signal "sn". The WVD coeffitions are stored
 * in "coefs", a complex matrix. The column represents time, and row
 * represents frequency.
 */
template <typename Type>
Matrix<Type> wvd( const Vector<Type> &sn )
{
    int N = sn.size(),
        dN = 2*N;

    Vector<Type> xn = fftInterp( sn, 2 );
    Vector<Type> fn( 3*dN );
    for( int i=dN; i<2*dN; ++i )
        fn[i] = xn[i-dN];

    Vector<Type> yn( dN );
    Matrix<Type> coefs( N, N );

    for( int n=1; n<=N; ++n )
    {
        for( int i=0; i<N; ++i )
            yn[i] = fn(dN+2*n+i) * fn(dN+2*n-i);
        for( int i=-N; i<0; ++i )
            yn[dN+i] = fn(dN+2*n+i) * fn(dN+2*n-i);

        coefs.setColumn( dyadDown(real(fft(yn)),0), n-1 );
    }

    return coefs;
}


/**
 * Compute WVD of 1D complex signal "cn". The WVD coeffitions are stored
 * in "coefs", a complex matrix. The column represents time, and row
 * represents frequency.
 */
template <typename Type>
Matrix<Type> wvd( const Vector< complex<Type> > &cn )
{
    int N = cn.size(),
        dN = 2*N;

    Vector< complex<Type> > xn = fftInterp( cn, 2 );
    Vector< complex<Type> > fn(3*dN);
    for( int i=dN; i<2*dN; ++i )
        fn[i] = xn[i-dN];

    Vector< complex<Type> > yn( dN );
    Matrix<Type> coefs( N, N );

    for( int n=1; n<=N; ++n )
    {
        for( int i=0; i<N; ++i )
            yn[i] = fn(dN+2*n+i) * conj(fn(dN+2*n-i));
        for( int i=-N; i<0; ++i )
            yn[dN+i] = fn(dN+2*n+i) * conj(fn(dN+2*n-i));

        coefs.setColumn( dyadDown(real(fft(yn)),0), n-1 );
    }

    return coefs;
}
