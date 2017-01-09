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
 *                               dwt-impl.h
 *
 * Implementation for DWT class.
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * constructors and destructor
 */
template<typename Type>
DWT<Type>::DWT( const string &wname ) : waveType(wname)
{
    if( waveType != "db4" )
    {
        cerr << "No such wavelet type!" << endl;
        exit(1);
    }

    getFilter(wname);
}

template<typename Type>
DWT<Type>::~DWT()
{
}


/**
 * get coefficients of filter bank
 */
template<typename Type>
inline void DWT<Type>::getFilter( const string &wname )
{
    if( wname == "db4" )
        db4Coefs( ld, hd, lr, hr );
}


/**
 * Initializing the length of coefficients.
 */
template <typename Type>
void DWT<Type>::lengthInit( int sigLength, int J )
{
    lenInfo.resize(J+2);

    int na = sigLength,
        nd = 0,
        total = 0;

    lenInfo[0] = sigLength;
    for( int j=1; j<=J; ++j )
    {
        // detial's length at level j
        nd = (na+hd.dim()-1)/2;

        // approx's length at level j
        na = (na+ld.dim()-1)/2;

        lenInfo[j] = nd;
        total += nd;
    }
    total += na;

    lenInfo[J+1] = na;
}


/**
 * Get the jth level approximation coefficients.
 */
template <typename Type>
Vector<Type> DWT<Type>::getApprox( const Vector<Type> &coefs )
{
    int J = lenInfo.dim()-2;
    Vector<Type> approx(lenInfo[J+1]);

    int start = 0;
    for( int i=1; i<=J; ++i )
        start += lenInfo[i];

    for( int i=0; i<lenInfo[J+1]; ++i )
        approx[i] = coefs[i+start];

    return approx;
}


/**
 * Get the Jth level detial coefficients.
 */
template <typename Type>
Vector<Type> DWT<Type>::getDetial( const Vector<Type> &coefs, int j )
{
    if( (j < 1) || (j > lenInfo.dim()-2) )
    {
        cerr << "Invalid level for getting detial coefficients!" << endl;
        return Vector<Type>(0);
    }

    Vector<Type> detial(lenInfo[j]);

    int start = 0;
    for( int i=1; i<j; ++i )
        start += lenInfo[i];

    for( int i=0; i<lenInfo[j]; ++i )
        detial[i] = coefs[i+start];

    return detial;
}


/**
 * Set the Jth level approximation coefficients.
 */
template <typename Type>
void DWT<Type>::setApprox( const Vector<Type> &approx, Vector<Type> &coefs )
{
    int J = lenInfo.dim()-2;
    if( approx.dim() != lenInfo[J+1] )
    {
        cerr << "Invalid length for setting approximation coefficient!"
             << endl;
        return;
    }

    // the approximation's start position in the coe
    int start = 0;
    for( int i=1; i<=J; ++i )
        start += lenInfo[i];

    for( int i=0; i<lenInfo[J+1]; ++i )
        coefs[i+start] = approx[i];
}


/**
 * Set the jth level's detial coefficients.
 */
template <typename Type>
void DWT<Type>::setDetial( const Vector<Type> &detial, Vector<Type> &coefs,
                           int j )
{
    if( detial.dim() != lenInfo[j] )
    {
        cerr << "Invalid length for setting detial coefficients!" << endl;
        return;
    }

    // the jth level detial's start position in the coe
    int start = 0;
    for( int i=1; i<j; ++i )
        start += lenInfo[i];

    for( int i=0; i<lenInfo[j]; ++i )
        coefs[i+start] = detial[i];
}


/**
 * Doing J levels decompose for input signal "signal".
 */
template <typename Type>
Vector<Type> DWT<Type>::dwt( const Vector<Type> &signal, int J )
{
    if( (lenInfo.size() != J+2) || (lenInfo[0] != signal.size()) )
    {
        cout << "first" << endl;
        lengthInit( signal.size(), J );
    }
//        lengthInit( signal.size(), J );

    int total = 0;
    for( int i=1; i<=J+1; ++i )
        total += lenInfo[i];

    Vector<Type> coefs(total);
    Vector<Type> approx = signal;
    Vector<Type> tmp;

    for( int j=1; j<=J; ++j )
    {
        // high frequency filtering
        tmp = conv( approx, hd );

        // stored the downsampled signal
        setDetial( dyadDown(tmp,1), coefs, j );

        // low frequency filtering
        tmp = conv( approx, ld );
        approx = dyadDown( tmp, 1 );
    }
    setApprox( approx, coefs );

    return coefs;
}


/**
 * Recover the jth level approximation signal.
 * The default parameter is "j=0", which means recover the original signal.
 */
template <typename Type>
Vector<Type> DWT<Type>::idwt( const Vector<Type> &coefs, int j )
{
    if( coefs.size() != (sum(lenInfo)-lenInfo[0]) )
    {
        cerr << "Invalid wavelet coeffcients!" << endl;
        return Vector<Type>(0);
    }

    int J = lenInfo.dim() - 2;
    if( (j < 0) || (j > J) )
    {
        cerr << "Invalid level for reconstructing signal!" << endl;
        return Vector<Type>(0);
    }

    Vector<Type> signal = getApprox( coefs );
    Vector<Type> detial;

    for( int i=J; i>j; --i )
    {
        detial = getDetial( coefs, i );

        // upsampling
        signal = dyadUp( signal, 0 );
        detial = dyadUp( detial, 0 );

        // recover the jth approximation
        signal = conv( signal, lr ) + conv( detial, hr );

        // cut off
        signal = wkeep( signal,lenInfo[i-1],"center" );
    }

    return signal;
}
