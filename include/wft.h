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
 *                                    wft.h
 *
 * Windowed Fourier Transform.
 *
 * These routines are designed for calculating WFT and IWFT of 1D signals.
 * The windows for forward and backword transform should be the same.
 *
 * In order to eliminate border effect, the input signal is extended by
 * three forms: zeros padded("zpd"), periodized extension("ppd") and
 * symetric extension("sym").
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef WFT_H
#define WFT_H


#include <string>
#include <complex>
#include <fft.h>
#include <matrix.h>
#include <utilities.h>


namespace splab
{

    template<typename Type>
    Matrix< complex<Type> > wft( const Vector<Type>&, const Vector<Type>&,
                                 const string &mode = "zpd" );

    template<typename Type>
    Vector<Type> iwft( const Matrix< complex<Type> >&, const Vector<Type>& );


    #include <wft-impl.h>

}
// namespace splab


#endif
// WFT_H
