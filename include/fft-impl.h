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
 *                                 fft-impl.h
 *
 * Implementation for FFT and IFFT interface.
 *
 * Zhang Ming, 2010-09, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Forward FFT algorithm.
 */
template<typename Type>
inline Vector< complex<Type> > fft( const Vector<Type> &xn )
{
    return fftr2c( xn );
}

template<typename Type>
inline Vector< complex<Type> > fft( const Vector< complex<Type> > &xn )
{
    return fftc2c( xn );
}


/**
 * Inverse FFT algorithm.
 */
template<typename Type>
inline Vector< complex<Type> > ifft( const Vector< complex<Type> > &Xk )
{
    return ifftc2c( Xk );
}


/**
 * Real to complex DFT of 1D signal.
 */
template<typename Type>
inline Vector< complex<Type> > fftr2c( const Vector<Type> &xn )
{
    Vector< complex<Type> > Xk( xn.size() );

    if( isPower2(xn.size())  )
    {
        FFTMR<Type> dft;
        dft.fft( xn, Xk );
    }
    else
    {
        FFTPF<Type> dft;
        dft.fft( xn, Xk );
    }

    return Xk;
}


/**
 * Complex to complex DFT of 1D signal.
 */
template<typename Type>
inline Vector< complex<Type> > fftc2c( const Vector< complex<Type> > &xn )
{
    Vector< complex<Type> > Xk( xn.size() );

    if( isPower2(xn.size())  )
    {
        for( int i=0; i<xn.size(); ++i )
            Xk[i] = xn[i];
        FFTMR<Type> dft;
        dft.fft( Xk );
    }
    else
    {
        FFTPF<Type> dft;
        dft.fft( xn, Xk );
    }

    return Xk;
}


/**
 * Complex to real IDFT of 1D signal.
 */
template<typename Type>
inline Vector<Type> ifftc2r( const Vector< complex<Type> > &Xk )
{
    Vector<Type> xn( Xk.size() );

    if( isPower2(xn.size())  )
    {
        Vector< complex<Type> > tmp( Xk );
        FFTMR<Type> dft;
        dft.ifft( tmp, xn );
    }
    else
    {
        FFTPF<Type> dft;
        dft.ifft( Xk, xn );
    }

    return xn;
}


/**
 * Complex to complex IDFT of 1D signal.
 */
template<typename Type>
inline Vector< complex<Type> > ifftc2c( const Vector< complex<Type> > &Xk )
{
    Vector< complex<Type> > xn( Xk.size() );

    if( isPower2(xn.size())  )
    {
        for( int i=0; i<xn.size(); ++i )
            xn[i] = Xk[i];
        FFTMR<Type> dft;
        dft.ifft( xn );
    }
    else
    {
        FFTPF<Type> dft;
        dft.ifft( Xk, xn );
    }

    return xn;
}
