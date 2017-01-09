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
 *                                 window.h
 *
 * Implementation for window function.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


/**
 * Get the specified window.
 */
template <typename Type>
Vector<Type> window( const string &wnName, int N, Type amp )
{
    if( wnName == "Rectangle" )
        return rectangle( N, amp );
    else if( wnName == "Bartlett" )
        return bartlett( N, amp );
    else if( wnName == "Hanning" )
        return hanning( N, amp );
    else if( wnName == "Hamming" )
        return hamming( N, amp );
    else if( wnName == "Blackman" )
        return blackman( N, amp );
    else
    {
        cerr << "No such type window!" << endl;
        return Vector<Type> (0);
    }
}

template <typename Type>
Vector<Type> window( const string &wnName, int N, Type alpha, Type amp )
{
    if( wnName=="Kaiser" )
        return kaiser( N, alpha, amp );
    else if( wnName=="Gauss" )
        return gauss( N, alpha, amp );
    else
    {
        cerr << "No such type window!" << endl;
        return Vector<Type> (0);
    }
}


/**
 * Calculates rectangle window coefficients.
 */
template <typename Type>
Vector<Type> rectangle( int N, Type amp )
{
    Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        win[i] = amp;
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates bartlett window coefficients.
 */
template <typename Type>
Vector<Type> bartlett( int N, Type amp )
{
    Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        win[i] = amp*2*i / (N-1);
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates hanning window coefficients.
 */
template <typename Type>
Vector<Type> hanning( int N, Type amp )
{
   Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        win[i] = amp * Type( 0.5 - 0.5*cos(TWOPI*i/(N-1)) );
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates hamming window coefficients.
 */
template <typename Type>
Vector<Type> hamming( int N, Type amp )
{
    Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        win[i] = amp * Type( 0.54 - 0.46*cos(TWOPI*i/(N-1.0)) );
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates hamming window coefficients.
 */
template <typename Type>
Vector<Type> blackman( int N, Type amp )
{
    Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        win[i] = amp * Type ( 0.42 - 0.50*cos(TWOPI*i/(N-1.0))
                              + 0.08*cos(2*TWOPI*i/(N-1.0)) );
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates hamming window coefficients.
 */
template <typename Type>
Vector<Type> kaiser( int N, Type alpha, Type amp )
{
    Vector<Type> win(N);

    for( int i=0; i<(N+1)/2; ++i )
    {
        Type beta = 2*alpha * Type( sqrt(i*(N-i-1.0))/(N-1.0) );
        win[i] = amp * I0(beta) / I0(alpha);
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * Calculates gauss window coefficients. "Alpha: is a optional parameter,
 * the default value is 2.5.
 */
template <typename Type>
Vector<Type> gauss( int N, Type alpha, Type amp )
{
    Vector<Type> win(N);
    Type center = (N-1)/Type(2);

    for( int i=0; i<(N+1)/2; ++i )
    {
        Type tmp = alpha*(i-center) / center;
        win[i] = amp * Type( exp(-0.5*tmp*tmp ) );
        win[N-1-i] = win[i];
    }

    return win;
}


/**
 * The zeroth N modified Bessel function of the first kind.
 */
template <typename Type>
Type I0( Type alpha )
{
    double  J = 1.0,
            K = alpha / 2.0,
            iOld = 1.0,
            iNew;
    bool    converge = false;

    // Use series expansion definition of Bessel.
    for( int i=1; i<MAXTERM; ++i )
    {
        J *= K/i;
        iNew = iOld + J*J;

        if( (iNew-iOld) < EPS )
        {
            converge = true;
            break;
        }
        iOld = iNew;
    }

    if( !converge )
        return Type(0);

    return Type(iNew);
}
