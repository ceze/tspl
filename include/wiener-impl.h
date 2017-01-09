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
 *                              wiener-impl.h
 *
 * Implementation for Wiener Filter.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * By given the observed signal "xn" and desired signal "dn", this routine
 * return the coefficients for Wiener filter with "p" order.
 */
template <typename Type>
Vector<Type> wienerFilter( const Vector<Type> &xn,
                           const Vector<Type> &dn, int p )
{
    int N = xn.size();
    assert( dn.size() == N );

    // auto-correlation and corss-correlation
    Vector<Type> Rxx = fastCorr( xn, "unbiased" );
    Vector<Type> Rdx = fastCorr( dn, xn, "unbiased" );

    Vector<Type> tn(p+1), bn(p+1);
    for( int i=0; i<=p; ++i )
    {
        tn[i] = Rxx[N-1+i];
        bn[i] = Rdx[N-1+i];
    }

    // solving Wiener-Hopf equations
    return levinson( tn, bn );
}


/**
 * One step Wiener predictor by using a "p" order filter. The input
 * signal "xn" should be much longer the "p" in order to have a more
 * precision estimation of the Correlation Matrix of "xn".
 */
template <typename Type>
Vector<Type> wienerPredictor( const Vector<Type> &xn, int p )
{
    int N = xn.size();

    // auto-correlation and corss-correlation
    Vector<Type> Rxx = fastCorr( xn, "unbiased" );
    Vector<Type> tn(p+1), bn(p+1), predictor(p);

    for( int i=0; i<=p; ++i )
        tn[i] = Rxx[N-1+i];
    bn(1) = Type(1.0);

    // solving Yule-Walker equations
    Vector<Type> tmp = levinson( tn, bn );

    for( int i=1; i<=p; ++i )
        predictor(i) = -tmp(i+1) / tmp(1);

    return predictor;
}
