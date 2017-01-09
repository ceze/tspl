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
 *                                  utilities.h
 *
 * Some usable routines converted from "Matlab", which are used in wavelet
 * transform and time-frequency analysis, such as "filp"(same to reverse),
 * "shift", "circshift", "fftshift", "dyadup", "wkeep", "wextend" and so on.
 *
 * Zhang Ming, 2010-01, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef UTILITIES_H
#define UTILITIES_H


#include <string>
#include <vector.h>
#include <fft.h>


namespace splab
{

    int mod( int, int );
    int ceil( int, int );

	template<typename Type> Vector<Type> reverse( const Vector<Type>& );
    template<typename Type> Vector<Type> flip( const Vector<Type>& );
    template<typename Type> Vector<Type> shift( const Vector<Type>& );
    template<typename Type> Vector<Type> cirshift( const Vector<Type>& );
    template<typename Type> Vector<Type> fftshift( const Vector<Type>& );

    template<typename Type> Vector<Type> dyadUp( const Vector<Type>&, int );
    template<typename Type> Vector<Type> dyadDown( const Vector<Type>&, int );
    template<typename Type> Vector<Type> fftInterp( const Vector<Type>&, int );
    template<typename Type>
    Vector< complex<Type> > fftInterp( const Vector< complex<Type> >&, int );

    template<typename Type>
    Vector<Type> wkeep( const Vector<Type>&, int, int );
    template<typename Type>
    Vector<Type> wkeep( const Vector<Type>&, int,
                        const string &direction="center" );
    template<typename Type>
    Vector<Type> wextend( const Vector<Type>&, int,
                          const string &direction="both",
                          const string &mode="zpd" );


    #include <utilities-impl.h>

}
// namespace splab


#endif
// UTILITIES_H
