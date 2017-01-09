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
 *                                    dgt.h
 *
 * Discrete Gabor Transform.
 *
 * These routines are designed for calculating discrete Gabor transform and
 * its inversion of 1D signals. In order to eliminate the border effect, the
 * input signal("signal") is extended by three forms: zeros padded("zpd"),
 * periodized extension("ppd") and symetric extension("sym").
 *
 * The analysis/synthesis function is given by users, and it's daul
 * (synthesis/analysis) function can be computed by "daul" routine. The over
 * sampling rate is equal to N/dM, where N denotes frequency sampling numbers
 * and dM denotes the time sampling interval.
 *
 * N and dM should can be devided evenly by the window length "Lw". The
 * recovered signal just has the elements from 1 to dM*floor(Ls/dM) of the
 * original signal. So you'd better let dM can be deviede evenly by  the
 * original signal length "Ls".
 *
 * Zhang Ming, 2010-03, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef DGT_H
#define DGT_H


#include <string>
#include <fft.h>
#include <matrix.h>
#include <linequs1.h>
//#include <linequs3.h>
#include <utilities.h>


namespace splab
{


    template<typename Type>
    Vector<Type> daul( const Vector<Type>&, int, int );

    template<typename Type>
    Matrix< complex<Type> > dgt( const Vector<Type>&,
                                 const Vector<Type>&,
                                 int, int,
                                 const string &mode = "zpd" );
    template<typename Type>
    Vector<Type> idgt( const Matrix< complex<Type> >&,
                       const Vector<Type>&,
                       int, int );


    #include <dgt-impl.h>

}
// namespace splab


#endif
// DGT_H
