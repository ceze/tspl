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
 *                             convolution-impl.h
 *
 * Implementation for linear convolution.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * convolution and ploynonal multiplication.
 */
template <typename Type>
Vector<Type> conv( const Vector<Type> &signal, const Vector<Type> &filter )
{
    if( signal.dim() < filter.dim() )
        return convolution( filter, signal );
    else
        return convolution( signal, filter );
}

template <typename Type>
Vector<Type> convolution( const Vector<Type> &signal, const Vector<Type> &filter )
{
    int sigLength = signal.dim();
    int filLength = filter.dim();
    assert( sigLength >= filLength );

    int length = sigLength + filLength - 1;
    Vector<Type> x(length);

    for( int i=1; i<=length; ++i )
    {
        x(i) = 0;
        if( i < filLength )
            for( int j=1; j<=i; ++j )
                x(i) += filter(j) * signal(i-j+1);
        else if( i <= sigLength )
            for( int j=1; j<=filLength; ++j )
                x(i) += filter(j) * signal(i-j+1);
        else
            for( int j=i-sigLength+1; j<=filLength; ++j )
                x(i) += filter(j) * signal(i-j+1);
    }
    return x;
}


/**
 * Fast convolution by FFT.
 */
template<typename Type>
Vector<Type> fastConv( const Vector<Type> &xn, const Vector<Type> &yn )
{
    int M = xn.dim(),
        N = yn.dim();

    Vector<Type> xnPadded = wextend( xn, N-1, "right", "zpd" ),
                 ynPadded = wextend( yn, M-1, "right", "zpd" );
    return ifftc2r( fft(xnPadded) * fft(ynPadded) );

//    Vector< complex<Type> > Zk = fft(xnPadded) * fft(ynPadded);
//    return ifftc2r(Zk);

//    return ifftc2r( fft(wextend(xn,N-1,"right","zpd")) * fft(wextend(yn,M-1,"right","zpd")) );
}
