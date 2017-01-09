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
 *                                 correlation.h
 *
 * The routines in this file estimate the cross-correlation sequence of a
 * random process. Autocorrelation is handled as a special case.
 *
 * c = corr(x,y,opt) returns the cross-correlation sequence in a length
 * 2*N-1 vector, where x and y are length N vectors (N>1). If x and y are
 * not the same length, the shorter vector is zero-padded to the length of
 * the longer vector.
 *
 * The parameter "opt" specifies a normalization option for the cross-
 * correlation, where 'opt' is
 * "none"       : to use the raw, unscaled cross-correlations (default);
 * "biased"     : biased estimate of the cross-correlation function;
 * "unbiased"   : unbiased estimate of the cross-correlation function.
 *
 * We use FFT computing the auto-corelation and cross-corelation functions
 * based on fallowing facts: for real functions,
 * R1[x(t),y(t)] = sum{ x(u)*y(u-t) } = Conv[x(t),y(-t)]
 * R2[x(t),y(t)] = sum{ x(u)*y(u+t) } = Conv[x(-t),y(t)]
 * And here we use the first defination.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef CORRELATION_H
#define CORRELATION_H


#include <convolution.h>
#include <utilities.h>


namespace splab
{

    template<typename Type> Vector<Type> corr( const Vector<Type>&,
                                               const string &opt="none" );
    template<typename Type> Vector<Type> corr( const Vector<Type>&,
                                               const Vector<Type>&,
                                               const string &opt="none" );

    template<typename Type> Vector<Type> fastCorr( const Vector<Type>&,
                                                   const string &opt="none" );
    template<typename Type> Vector<Type> fastCorr( const Vector<Type>&,
                                                   const Vector<Type>&,
                                                   const string &opt="none" );

    template<typename Type> static void biasedProcessing( Vector<Type> &,
                                                          const string &opt );


    #include <correlation-impl.h>

}
// namespace splab


#endif
// CORRELATION_H
