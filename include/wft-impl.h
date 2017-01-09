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
 *                                wft-impl.h
 *
 * Implementation for Windowed Fourier Transform.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Compute WFT of 1D signal("xn") using "wn" as window. The time-frequency
 * coeffitions are stored in "coefs". The column represents time, and row
 * represents frequency.
 */
template <typename Type>
Matrix< complex<Type> > wft( const Vector<Type> &xn, const Vector<Type> &wn,
                             const string &mode )
{
    int Lx = xn.size(),
        Lw = wn.size();

    // extends the input signal
    Vector<Type> tmp = wextend( xn, Lw/2, "both", mode );

    Matrix< complex<Type> > coefs( Lw, Lx );
    Vector<Type> sn( Lw );
    Vector< complex<Type> > Sk( Lw );

    for( int i=0; i<Lx; ++i )
    {
        // intercept xn by wn function
        for( int j=0; j<Lw; ++j )
            sn[j] = tmp[i+j] * wn[j];
        Sk = fft(sn);

        // compute the Foureier transform
        coefs.setColumn( Sk, i );
    }

    return coefs;
}


/**
 * Compute the inverse windowed Fourier transform of "coefs". The window
 * "wn" should be the same as forward transform. The reconstruction signal
 * is stored in "xn".
 */
template <typename Type>
Vector<Type> iwft( const Matrix< complex<Type> > &coefs,
                   const Vector<Type> &wn )
{
    int Lw = wn.size(),
        Lx = coefs.cols();

    Vector<Type> xn(Lx);
    Matrix<Type> tmp( Lw, Lx );
    Vector< complex<Type> > Sk( Lw );
    Vector<Type> sn( Lw );

    // compute the inverse Fourier transform of coefs
    for( int i=0; i<Lx; ++i )
    {
        Sk = coefs.getColumn(i);
        sn = ifftc2r(Sk);
        tmp.setColumn( sn, i );
    }

    int mid = Lw / 2;
    for( int i=0; i<Lx; ++i )
        xn[i] = tmp[mid][i] / wn[mid];

    return xn;
}
