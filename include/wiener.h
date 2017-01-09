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
 *                                   wiener.h
 *
 * Wiener Filter.
 *
 * The goal of the Wiener filter is to filter out noise that has corrupted
 * a signal. It is based on a statistical approach.
 *
 * Wiener filters are characterized by the following:
 * Assumption:  signal and (additive) noise are stationary linear stochastic
 *              processes with known spectral characteristics or known
 *              autocorrelation and cross-correlation.
 * Requirement: the filter must be physically realizable/causal (this
 *              requirement can be dropped, resulting in a non-causal solution).
 * Performance: criterion: minimum mean-square error (MMSE).
 *
 * And The correlation matrix is a symmetric Toeplitz matrix, so we can use the
 * efficient algorithm, Levinson-Durbin algorithm, to solve the Wiener-Hopf
 * Equations.
 *
 * Zhang Ming, 2010-10, Xi'an Jiaotong University.
 *****************************************************************************/


#ifndef WIENER_H
#define WIENER_H


#include <vector.h>
#include <correlation.h>
#include <levinson.h>


namespace splab
{

    template<typename Type>
    Vector<Type> wienerFilter( const Vector<Type>&,
                               const Vector<Type>&, int );

    template<typename Type>
    Vector<Type> wienerPredictor( const Vector<Type>&, int );


    #include <wiener-impl.h>

}
// namespace splab


#endif
// WIENER_H
